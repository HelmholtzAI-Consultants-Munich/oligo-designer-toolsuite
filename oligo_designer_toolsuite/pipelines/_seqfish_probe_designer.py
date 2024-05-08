import os
import sys

sys.path.append(os.path.dirname(os.getcwd()))
import logging
import warnings

import yaml

from ..database import (
    CustomGenomicRegionGenerator,
    EnsemblGenomicRegionGenerator,
    NcbiGenomicRegionGenerator,
    OligoDatabase,
    ReferenceDatabase,
)
from ..oligo_efficiency_filter import AverageSetScoring, SeqFISHOligoScoring
from ..oligo_property_filter import (
    ConsecutiveRepeats,
    GCContent,
    MaskedSequences,
    PropertyFilter,
)
from ..oligo_selection import OligosetGenerator, padlock_heuristic_selection
from ..oligo_specificity_filter import Blastn, ExactMatches, SpecificityFilter

# from ._seqfish_readout_probe_designer import SeqFISHReadoutProbeDesigner
# from ..sequence_design._barcoding_creation import BarcodingCreator
# from ..sequence_design._seqFISH_sequence import SeqfishProbesCreator


class SeqFISHProbeDesigner:
    """
    This class is designed to generate primary probes of SeqFISH+ experiment
    """

    def __init__(self, config_file) -> None:
        """
        Constructor method
        :param config_file: path to yaml configuration file
        :type config_file: str
        """
        self.config_file = config_file
        with open(config_file, "r") as yaml_file:
            self.config = yaml.safe_load(yaml_file)
        self.dir_output = os.path.join(os.path.dirname(os.getcwd()), self.config["dir_output"])

    def create_probes(self):
        """ "
        Method to create probes using config file.
        """
        logging.info(f"Create transcriptome " + str(self.config["source"]))

        if self.config["source"] == "custom":
            self.region_generator = CustomGenomicRegionGenerator(
                annotation_file=self.config["file_annotation"],
                sequence_file=self.config["file_sequence"],
                files_source=self.config["files_source"],
                species=self.config["species"],
                annotation_release=self.config["annotation_release"],
                genome_assembly=self.config["genome_assembly"],
                dir_output=self.dir_output,
            )
        elif self.config["source"] == "ncbi":
            self.region_generator = NcbiGenomicRegionGenerator(
                taxon=self.config["taxon"],
                species=self.config["species"],
                annotation_release=self.config["annotation_release"],
                dir_output=self.dir_output,
            )
        elif self.config["source"] == "ensembl":
            self.region_generator = EnsemblGenomicRegionGenerator(
                species=self.config["species"],
                annotation_release=self.config["annotation_release"],
                dir_output=self.dir_output,
            )
        logging.info(f"Create file_transcriptome")
        self.file_transcriptome = self.region_generator.generate_transcript_reduced_representation(
            include_exon_junctions=True,
            exon_junction_size=2 * self.config["oligo_length_max"],
        )
        logging.info(f"Initialisation is finished")
        # CREATE OLIGO DATABASE
        self.oligo_database = OligoDatabase(
            file_fasta=self.file_transcriptome,
            oligo_length_min=self.config["oligo_length_min"],
            oligo_length_max=self.config["oligo_length_max"],
            min_oligos_per_region=self.config["min_oligos_per_gene"],
            files_source=self.region_generator.files_source,
            species=self.region_generator.species,
            annotation_release=self.region_generator.annotation_release,
            genome_assembly=self.region_generator.genome_assembly,
            n_jobs=2,
            dir_output=self.dir_output,
        )

        if self.config["file_genes"] is None:
            warnings.warn(
                "No file containing the genes was provided, all the genes are ussed to generate the probes. This chioce can use a lot of resources."
            )
            genes = None
        else:
            with open(self.config["file_genes"]) as handle:
                lines = handle.readlines()
                genes = [line.rstrip() for line in lines]

        # generate the oligo sequences from gene transcripts
        self.oligo_database.create_database(region_ids=genes)
        logging.info(f"Oligo database created")

        # WRITE INTERMIDIATE RESULTS
        if self.config["write_intermediate_steps"]:
            file_database = self.oligo_database.write_database(filename="oligo_database_initial.txt")

        # PROPERTY FILTERS
        # Filters: MaskedSequences, GCContent, Prohibited Sequences
        masked_sequences = MaskedSequences()
        gc_content = GCContent(
            GC_content_min=self.config["GC_content_min"],
            GC_content_max=self.config["GC_content_max"],
        )
        proh_seq = ConsecutiveRepeats(num_consecutive=self.config["number_consecutive"])
        filters = [masked_sequences, proh_seq, gc_content]
        property_filter = PropertyFilter(
            filters=filters,
            write_regions_with_insufficient_oligos=self.config["write_removed_genes"],
        )
        self.oligo_database = property_filter.apply(
            oligo_database=self.oligo_database, n_jobs=self.config["n_jobs"]
        )
        self.oligo_database.remove_regions_with_insufficient_oligos(
            pipeline_step="after applying property filters"
        )
        if self.config["write_intermediate_steps"]:
            file_database = self.oligo_database.write_database(filename="oligo_database_property_filter.txt")

        logging.info(f"Property filters applied")

        # SPECIFICITY FILTERS
        # Filters: ExactMatches and Blastn

        self.dir_specificity = os.path.join(
            self.dir_output, "specificity_temporary"
        )  # folder where the temporary files will be written
        self.reference = ReferenceDatabase(
            file_fasta=self.file_transcriptome,
            files_source=self.config["files_source"],
            species=self.config["species"],
            annotation_release=self.config["annotation_release"],
            genome_assembly=self.config["genome_assembly"],
            dir_output=self.dir_output,
        )
        self.reference.load_fasta_into_database()
        logging.info(f"Reference DB created")
        exact_mathces = ExactMatches(dir_specificity=self.dir_specificity)
        blastn = Blastn(
            dir_specificity=self.dir_specificity,
            word_size=self.config["word_size"],
            percent_identity=self.config["percent_identity"],
            coverage=self.config["coverage"],
            strand=self.config["strand"],
            # strand='plus',
        )
        filters = [exact_mathces, blastn]
        specificity_filter = SpecificityFilter(
            filters=filters,
            write_regions_with_insufficient_oligos=self.config["write_removed_genes"],
        )
        self.oligo_database = specificity_filter.apply(
            oligo_database=self.oligo_database,
            reference_database=self.reference,
            n_jobs=self.config["n_jobs"],
        )
        if self.config["write_intermediate_steps"]:
            file_database = self.oligo_database.write_database(
                filename="oligo_database_specificity_filter.txt"
            )

        logging.info(f"Specificity filters applied")

        # READOUT PROBES GENERATOR
        # Readout probes are created, returned as a dictionary self.readout_probes and also stored in the file
        # inside self.dir_output
        readouts_generator = SeqFISHReadoutProbeDesigner(self.config, blastn, self.reference, self.dir_output)
        self.readout_probes = readouts_generator.create_readout_probes()
        logging.info(f"Readout probes are created")

        # BARCODES CREATOR
        # Barcodes for each gene are being created here, returned as a dictionary and
        # also stored in a .txt file inside self.dir_output
        barcodes_creator = BarcodingCreator(
            self.config["num_pseudocolors"], list(self.oligo_database.database.keys())
        )
        barcodes_for_genes = barcodes_creator.create_barcodes()
        output_file_barcodes = os.path.join(self.dir_output, "barcodes_for_each_gene.txt")
        f = open(output_file_barcodes, "w")
        for i in barcodes_for_genes.keys():
            S = "["
            for j in barcodes_for_genes[i]:
                S = S + str(j) + ", "
            S = S + "]\n"
            f.write(i + " : " + S)
        logging.info(f"Barcodes are assigned")

        # Using barcodes total probes are being assembled here and stored in oligoDB
        probes_creator = SeqfishProbesCreator()
        self.oligo_database.database = probes_creator.create_probes(
            self.oligo_database.database, self.readout_probes, barcodes_for_genes
        )
        logging.info(f"Total probes are built")

        # CROSS_HYBRIDIZATION CHECK
        # Check for cross-hybridization is performed here, ExactMatches and Blastn filters are being used
        # Blastn filter is used with strand "minus"
        self.oligo_database.write_fasta_from_database(filename="fasta_from_our_db")
        self.ref_db_cross_hybr = ReferenceDatabase(
            file_fasta=self.dir_output + "/oligo_database/fasta_from_our_db.fna"
        )
        exact_mathces = ExactMatches(dir_specificity=self.dir_specificity)
        blastn_cross_hybr = Blastn(
            dir_specificity=self.dir_specificity,
            word_size=self.config["word_size"],
            percent_identity=self.config["percent_identity"],
            coverage=self.config["coverage"],
            # THE MOST IMPORTANT PART HERE IS STRAND = "MINUS"
            strand="minus",
        )
        filters = [exact_mathces, blastn_cross_hybr]
        specificity_filter = SpecificityFilter(
            filters=filters,
            write_regions_with_insufficient_oligos=self.config["write_removed_genes"],
        )
        self.oligo_database = specificity_filter.apply(
            oligo_database=self.oligo_database,
            reference_database=self.ref_db_cross_hybr,
            n_jobs=self.config["n_jobs"],
        )

        logging.info(f"Cross-hybridization check is performed")

        # The best non-overlapping set of probes is built here
        # SeqFISH scoring function is used, see its descr in SeqFISHOligoScoring
        oligos_scoring = SeqFISHOligoScoring(
            GC_content_min=self.config["GC_content_min"],
            GC_content_opt=self.config["GC_content_opt"],
            GC_content_max=self.config["GC_content_max"],
            GC_weight=self.config["GC_weight"],
        )
        set_scoring = AverageSetScoring()
        oligoset_generator = OligosetGenerator(
            oligoset_size=self.config["oligoset_size"],
            min_oligoset_size=self.config["min_oligoset_size"],
            oligos_scoring=oligos_scoring,
            set_scoring=set_scoring,
            heurustic_selection=padlock_heuristic_selection,
            write_regions_with_insufficient_oligos=self.config["write_removed_genes"],
        )
        self.oligo_database = oligoset_generator.apply(
            oligo_database=self.oligo_database,
            n_sets=self.config["n_sets"],
            n_jobs=self.config["n_jobs"],
        )
        if self.config["write_intermediate_steps"]:
            self.oligo_database.write_oligosets(dir_oligosets="oligosets")
        logging.info(f"Oligoset generation is finished. Done!")
        return self.oligo_database
