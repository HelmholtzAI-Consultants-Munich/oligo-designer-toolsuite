import logging
import os
import shutil

# from typing_extensions import Literal # Python 3.7 or below
from typing import List, Literal, Tuple

from ..database import OligoDatabase, ReferenceDatabase
from ..oligo_efficiency_filter import AverageSetScoring, SeqFISHOligoScoring
from ..oligo_property_filter import (
    GCContentFilter,
    HardMaskedSequenceFilter,
    HomopolymericRunsFilter,
    PropertyFilter,
)
from ..oligo_selection import OligosetGenerator, padlock_heuristic_selection
from ..oligo_specificity_filter import Blastn, ExactMatches, SpecificityFilter
from ._base_oligo_designer import BaseOligoDesigner


from pathlib import Path

import yaml
from Bio.Seq import Seq
import numpy as np

from ..utils import create_seqfish_plus_barcodes

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


def create_seqfish_plus_barcodes(n_pseudocolors: int, seed: int, num_genes: int):
    """
    Function to create barcodes for each gene
    :return: dictionary of barcodes for each gene; { gene : barcode }
    :rtype: dict (str : list of 4 int)
    """
    output = dict()
    all_possible_barcodes = [
        [i, (i + j + k) % n_pseudocolors, j, k]
        for i in range(n_pseudocolors)
        for j in range(n_pseudocolors)
        for k in range(n_pseudocolors)
    ]
    l = len(all_possible_barcodes)
    arr = np.arange(0, l, 1)
    np.random.seed(seed)
    barcodes = np.random.choice(arr, size=num_genes, replace=False)
    output = [all_possible_barcodes[barcodes[i]] for i in range(num_genes)]
    return output


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


class SeqfishPlusProbeDesigner(BaseOligoDesigner):
    """_summary_

    Args:
        BaseOligoDesigner (_type_): _description_
    """

    # 1
    def filter_probes_by_property(
        self, oligo_database, GC_content_min, GC_content_max, number_consecutive, n_jobs
    ):
        masked_sequences = HardMaskedSequenceFilter()
        gc_content = GCContentFilter(GC_content_min=GC_content_min, GC_content_max=GC_content_max)
        proh_seq = HomopolymericRunsFilter(base=["A", "C", "T", "G"], n=number_consecutive)
        filters = [masked_sequences, proh_seq, gc_content]
        property_filter = PropertyFilter(filters=filters)
        oligo_database = property_filter.apply(oligo_database=oligo_database, n_jobs=n_jobs)
        if self.write_intermediate_steps:
            file_database = oligo_database.write_database(filename="oligo_database_property_filter.txt")
        return oligo_database, file_database

    # 2
    def filter_probes_by_specificity(
        self,
        oligo_database: OligoDatabase,
        probe_length_max: int,
        region_reference: Literal["cds", "reduced_representation", "genome"],
        word_size: int,
        percent_identity: float,
        coverage: float,
        strand: str,
        n_jobs: int,
    ) -> Tuple[OligoDatabase, str]:
        """_summary_

        Args:
            oligo_database (OligoDatabase): _description_
            probe_length_max (int): _description_
            word_size (int): _description_
            percent_identity (float): _description_
            coverage (float): _description_
            strand (str): _description_
            n_jobs (int): _description_

        Returns:
            Tuple[OligoDatabase, str]: _description_
        """

        num_genes_before, num_probes_before = self._get_probe_database_info(oligo_database.database)
        dir_specificity = os.path.join(
            self.dir_output, "specificity_temporary"
        )  # folder where the temporary files will be written

        if region_reference == "reduced_representation":
            file_transcriptome_ref = self.region_generator.generate_transcript_reduced_representation(
                include_exon_junctions=True, exon_junction_size=probe_length_max
            )
        elif region_reference == "genome":
            file_transcriptome_ref = self.region_generator.generate_genome()
        elif region_reference == "cds":
            file_transcriptome_ref = self.region_generator.generate_CDS_reduced_representation(
                include_exon_junctions=True, exon_junction_size=probe_length_max
            )
        self.reference = ReferenceDatabase(
            file_fasta=file_transcriptome_ref,
            metadata=self.metadata,
            dir_output=self.dir_output,
        )
        self.reference.load_fasta_into_database()
        logging.info(f"Reference DB created")
        exact_matches = ExactMatches(dir_specificity=dir_specificity)
        blastn = Blastn(
            dir_specificity=dir_specificity,
            word_size=word_size,  # 15 (config)
            percent_identity=percent_identity,  # 100 (config)
            coverage=coverage,  # alignment length instead (config)
            strand=strand,
        )
        filters = [exact_matches, blastn]
        specificity_filter = SpecificityFilter(filters=filters)

        oligo_database = specificity_filter.apply(
            oligo_database=oligo_database,
            reference_database=self.reference,
            n_jobs=n_jobs,
        )
        if self.write_intermediate_steps:
            file_database = oligo_database.write_database(filename="oligo_database_specificity_filter.txt")

        num_genes_after, num_probes_after = self._get_probe_database_info(oligo_database.database)

        logging.info(
            f"Step - Filter Probes by Specificity: the database contains {num_probes_after} probes from {num_genes_after} genes, while {num_probes_before - num_probes_after} probes and {num_genes_before - num_genes_after} genes have been deleted in this step."
        )

        shutil.rmtree(dir_specificity)
        return oligo_database, file_database

    # 5
    def design_readout_probes(
        self,
        GC_content_min: float,
        GC_content_max: float,
        number_consecutive: int,
        length: int,
        blast_filter: Blastn,
        sequence_alphabet: List[str],
        num_probes: int,
        reference_DB: ReferenceDatabase,
    ) -> OligoDatabase:
        """_summary_

        Args:
            GC_content_min (float): _description_
            GC_content_max (float): _description_
            number_consecutive (int): _description_
            length (int): _description_
            blast_filter (Blastn): _description_
            sequence_alphabet (List[str]): _description_
            num_probes (int): _description_
            reference_DB (ReferenceDatabase): _description_

        Returns:
            OligoDatabase: _description_
        """

        property_filters = [GCContentFilter(GC_content_min=GC_content_min, GC_content_max=GC_content_max)]

        # Reuse specificity filter and cross hybridization check (different parameters)
        specificity_filters = [blast_filter]

        readout_database = OligoDatabase(dir_output=self.dir_output)

        # TODO make sure to generate enough random probes
        # We multiply by 20 to make sure we have enough probes
        readout_database.create_random_database(length, num_probes * 20, sequence_alphabet=sequence_alphabet)

        property_filter = PropertyFilter(property_filters)
        readout_database = property_filter.apply(readout_database)

        specificity_filter = SpecificityFilter(specificity_filters)
        readout_database = specificity_filter.apply(readout_database, reference_DB)
        return readout_database

    # 3 (target)
    def check_cross_hybridization(
        self,
        oligo_database,
        dir_specificity,
        word_size,
        percent_identity,
        coverage,
        n_jobs,
    ):
        # CROSS_HYBRIDIZATION CHECK
        # Check for cross-hybridization is performed here, ExactMatches and Blastn filters are being used
        # Blastn filter is used with strand "minus"
        oligo_database.write_fasta_from_database(filename="fasta_from_our_db")
        ref_db_cross_hybr = ReferenceDatabase(
            file_fasta=self.dir_output + "/oligo_database/fasta_from_our_db.fna"
        )
        blastn_cross_hybr = Blastn(
            dir_specificity=dir_specificity,  # 15 (config) see specificity filter
            word_size=word_size,
            percent_identity=percent_identity,
            coverage=coverage,
            strand="minus",  # To be verified
        )
        filters = [blastn_cross_hybr]
        specificity_filter = SpecificityFilter(filters=filters)
        oligo_database = specificity_filter.apply(
            oligo_database=oligo_database,
            reference_database=ref_db_cross_hybr,
            n_jobs=n_jobs,
        )
        if self.write_intermediate_steps:
            file_database = oligo_database.write_database(filename="cross_hybridation.txt")

        logging.info(f"Cross-hybridization check is performed")
        return oligo_database, file_database

    # 4
    def get_oligosets(
        self,
        oligo_database,
        GC_content_opt,
        oligoset_size,
        min_oligoset_size,
        oligos_scoring,
        n_sets,
        n_jobs,
        max_oligos,
    ):
        oligos_scoring = SeqFISHOligoScoring(GC_content_opt=GC_content_opt)
        set_scoring = AverageSetScoring()

        # distance between oligos = 2 (config)
        oligoset_generator = OligosetGenerator(
            oligoset_size=oligoset_size,
            min_oligoset_size=min_oligoset_size,
            oligos_scoring=oligos_scoring,
            set_scoring=set_scoring,
            heurustic_selection=padlock_heuristic_selection,
            max_oligos=max_oligos,
        )
        oligo_database = oligoset_generator.apply(oligo_database=oligo_database, n_sets=n_sets, n_jobs=n_jobs)
        if self.write_intermediate_steps:
            dir_oligosets = oligo_database.write_oligosets(folder="oligosets")
            file_database = oligo_database.write_database(filename="oligo_database_oligosets.txt")
        logging.info(f"Oligoset generation is finished. Done!")

        return oligo_database, file_database, dir_oligosets


class SeqfishProbesCreator:
    """
    This class is used to assemble probes using primary probes, readout probes and barcodes, that were designed for each gene.
    """

    def __init__(self, dir_output: str):
        self.dir_output = dir_output

    def assemble_probes(
        self,
        oligo_database,
        readout_database,
        n_pseudocolors,
        seed,
    ):
        barcodes = create_seqfish_plus_barcodes(
            n_pseudocolors=n_pseudocolors,
            seed=seed,
            num_genes=len(oligo_database.database.keys()),
        )
        readout_sequences = readout_database.to_sequence_list()
        for i in oligo_database.keys():
            barcode = barcodes[i]
            left = readout_sequences[barcode[0]] + readout_sequences[barcode[1]]
            right = readout_sequences[barcode[2]] + readout_sequences[barcode[3]]

            for j in oligo_database[i].keys():
                seq = str(oligo_database[i][j]["sequence"])
                seq = left + " T " + seq  # T is a spacer, to be outputted separately
                seq = seq + " T " + right
                oligo_database[i][j]["sequence"] = Seq(seq)
        return oligo_database

    def write_final_probes(
        self,
        oligo_database,
        readout_database,
        n_pseudocolors,
        seed,
    ):
        barcodes = create_seqfish_plus_barcodes(
            n_pseudocolors=n_pseudocolors,
            seed=seed,
            num_genes=len(oligo_database.database.keys()),
        )
        readout_sequences = readout_database.to_sequence_list()

        regions = oligo_database.database.keys()
        yaml_dict = {}

        for region in regions:
            barcode = barcodes[region]
            left = readout_sequences[barcode[0]] + readout_sequences[barcode[1]]
            right = readout_sequences[barcode[2]] + readout_sequences[barcode[3]]

            yaml_dict[region] = {}
            for i, oligo in enumerate(oligo_database[region].keys()):
                yaml_dict[region][f"{region}_oligo_{i+1}"] = {}
                seq = str(oligo_database[region][oligo]["sequence"])
                yaml_dict[region][f"{region}_oligo_{i+1}"]["id"] = oligo
                yaml_dict[region][f"{region}_oligo_{i+1}"]["region"] = region
                yaml_dict[region][f"{region}_oligo_{i+1}"]["seqfish_plus_full_probe"] = (
                    left + "T" + seq + "T" + right
                )
                yaml_dict[region][f"{region}_oligo_{i+1}"].update(oligo_database[region][oligo])
                yaml_dict[region][f"{region}_oligo_{i+1}"]["start"] = ",".join(
                    f"{start}" for start in yaml_dict[region][f"{region}_oligo_{i+1}"].pop("start")
                )
                yaml_dict[region][f"{region}_oligo_{i+1}"]["end"] = ",".join(
                    f"{end}" for end in yaml_dict[region][f"{region}_oligo_{i+1}"].pop("end")
                )
                yaml_dict[region][f"{region}_oligo_{i+1}"]["target_sequence"] = yaml_dict[region][
                    f"{region}_oligo_{i+1}"
                ].pop("sequence")

                yaml_dict[region][f"{region}_oligo_{i+1}"].update(
                    {f"readout_sequece_{i+1}": readout_sequences[barcode[i]] for i in range(4)}
                )
        # save
        probes_dir = os.path.join(self.dir_output, "final_seqfish_plus_sequences")
        Path(probes_dir).mkdir(parents=True, exist_ok=True)
        with open(os.path.join(probes_dir, "seqfish_plus_probes.yml"), "w") as outfile:
            yaml.dump(yaml_dict, outfile, default_flow_style=False, sort_keys=False)

        # Create order file
        yaml_order = {}
        for region in yaml_dict:
            yaml_order[region] = {}
            for oligo_id in yaml_dict[region]:
                yaml_order[region][oligo_id] = {}
                yaml_order[region][oligo_id]["merfish_probe_full_sequence"] = yaml_dict[region][oligo_id][
                    "merfish_probe_full_sequence"
                ]
                yaml_order[region][oligo_id]["readout_probe_1"] = yaml_dict[region][oligo_id][
                    "readout_probe_1"
                ]
                yaml_order[region][oligo_id]["readout_probe_2"] = yaml_dict[region][oligo_id][
                    "readout_probe_2"
                ]
        with open(os.path.join(probes_dir, "merfish_probes_order.yml"), "w") as outfile:
            yaml.dump(yaml_order, outfile, default_flow_style=False, sort_keys=False)

        # Create readout probe file
        yaml_readout = {}
        yaml_readout["Bit"] = "Readout Probe"
        for i in range(num_bits):
            yaml_readout[str(i + 1)] = readout_probes[i] + "/3Cy5Sp"
        with open(os.path.join(probes_dir, "merfish_readout_probes.yml"), "w") as outfile:
            yaml.dump(yaml_readout, outfile, default_flow_style=False, sort_keys=False)

        # Create codebook file
        yaml_codebook = {}
        for gene_idx, gene in enumerate(genes):
            yaml_codebook[gene] = self.code[gene_idx]

        for i in range(self.n_blanks):
            yaml_codebook[f"blank_barcode_{i + 1}"] = self.code[n_genes + i]
        with open(os.path.join(probes_dir, "merfish_codebook.yml"), "w") as outfile:
            yaml.dump(yaml_codebook, outfile, default_flow_style=False, sort_keys=False)
