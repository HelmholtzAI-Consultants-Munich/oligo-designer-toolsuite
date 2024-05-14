import inspect
import logging
import os
import sys
from datetime import datetime
from pathlib import Path

# from typing_extensions import Literal # Python 3.7 or below
from typing import Literal

from oligo_designer_toolsuite.database import OligoDatabase
from oligo_designer_toolsuite.sequence_generator import (
    CustomGenomicRegionGenerator,
    EnsemblGenomicRegionGenerator,
    NcbiGenomicRegionGenerator,
    OligoSequenceGenerator,
)


class BaseOligoDesigner:
    """ """

    def __init__(
        self,
        dir_output: str = "output",
        log_name: str = "oligo_designer",
        write_removed_genes: bool = True,
        write_intermediate_steps: bool = True,
    ):
        """Constructor"""
        ##### store parameters #####
        self.dir_output = dir_output
        Path(self.dir_output).mkdir(parents=True, exist_ok=True)

        self.write_removed_genes = write_removed_genes
        self.write_intermediate_steps = write_intermediate_steps
        self.log_name = log_name

        ##### setup logger #####
        timestamp = datetime.now()
        file_logger = os.path.join(
            self.dir_output,
            f"log_{log_name}_{timestamp.year}-{timestamp.month}-{timestamp.day}-{timestamp.hour}-{timestamp.minute}.txt",
        )
        logging.getLogger("log_name")
        logging.basicConfig(
            format="%(asctime)s [%(levelname)s] %(message)s",
            level=logging.NOTSET,
            handlers=[logging.FileHandler(file_logger), logging.StreamHandler()],
        )
        logging.captureWarnings(True)

        ##### log parameters #####
        logging.info("Parameters Init:")
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        parameters = {i: values[i] for i in args}
        self._log_parameters(parameters)

        ##### initialize annotation parameters #####
        self.annotation_file = None
        self.sequence_file = None
        self.metadata = {}

    def _log_parameters(self, parameters):
        """Log function parameters.

        :param parameters: Dict with parameter name : parameter value pairs
        :type parameters: dict
        """
        for key, value in parameters.items():
            if key != "self":
                logging.info(f"{key} = {value}")

    def _get_oligo_database_info(self, oligo_database: dict):
        """Count the number of oligos and genes in the database.

        :param oligo_database: Database with oligos.
        :type oligo_database: dict
        :return: Number of genes and oligos in the database.
        :rtype: int, int
        """
        genes = oligo_database.keys()
        num_genes = len(genes)
        num_oligos = 0
        for gene in genes:
            num_oligos += len(oligo_database[gene].keys())

        return num_genes, num_oligos

    def _get_oligo_length_min_max_from_database(self, oligo_database: dict):
        """Get minimum and maximum length of oligos stored in the oligo database.

        :param oligo_database: Database with oligos.
        :type oligo_database: dict
        :return: Min and max length of oligos
        :rtype: int, int
        """
        oligo_length_min = sys.maxsize
        oligo_length_max = 0

        for region in oligo_database.keys():
            for oligo in oligo_database[region].keys():
                length = oligo_database[region][oligo]["length"]
                if length < oligo_length_min:
                    oligo_length_min = length
                if length > oligo_length_max:
                    oligo_length_max = length

        return oligo_length_min, oligo_length_max

    def load_annotations(
        self,
        source: str,
        source_params: dict,
    ):
        """Load annotations for specified source. Source can be either "ncbi", "ensemble" or "custom".

        If "ncbi" is choosen the source_params need to contain the keys:
        - "taxon": taxon of the species, valid taxa are: archaea, bacteria, fungi, invertebrate, mitochondrion, plant, plasmid, plastid, protozoa, vertebrate_mammalian, vertebrate_other, viral
        - "species": species name in NCBI download format, e.g. 'Homo_sapiens' for human; see [here](https://ftp.ncbi.nlm.nih.gov/genomes/refseq/) for available species name
        - "annotation_release": release number (e.g. 109 or 109.20211119 for ncbi) of annotation or 'current' to use most recent annotation release. Check out release numbers for NCBI at ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/annotation_releases/

        If "ensembl" is choosen the source_params need to contain the keys:
        - "species": species name in ensembl download format, e.g. 'homo_sapiens' for human; see http://ftp.ensembl.org/pub/release-108/gtf/ for available species names
        - "annotation_release": release number of annotation, e.g. 'release-108' or 'current' to use most recent annotation release. Check out release numbers for ensemble at ftp.ensembl.org/pub/

        If "custom" is choosen the source_params need to contain the keys:
        - "file_annotation": GTF file with gene annotation
        - "file_sequence": FASTA file with genome sequence
        - "files_source": original source of the genomic files -> optional, i.e. can be assigned None
        - "species": species of provided annotation, leave empty if unknown -> optional, i.e. can be assigned None
        - "annotation_release": release number of provided annotation, leave empty if unknown -> optional, i.e. can be assigned None
        - "genome_assembly": genome assembly of provided annotation, leave empty if unknown -> optional, i.e. can be assigned None

        :param source: Indicate from where the annotation files will be loaded. Options:  'ncbi', 'ensembl', 'custom'.
        :type source: str
        :param source_params: Parameters for loading annotations. See above for details.
        :type source_params: dict
        """
        ##### log parameters #####
        logging.info("Parameters Load Annotations:")
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        parameters = {i: values[i] for i in args}
        self._log_parameters(parameters)

        ##### loading annotations from different sources #####
        if source == "ncbi":
            # dowload the fasta files formthe NCBI server
            self.region_generator = NcbiGenomicRegionGenerator(
                taxon=source_params["taxon"],
                species=source_params["species"],
                annotation_release=source_params["annotation_release"],
                dir_output=self.dir_output,
            )
        elif source == "ensembl":
            # dowload the fasta files formthe NCBI server
            self.region_generator = EnsemblGenomicRegionGenerator(
                species=source_params["species"],
                annotation_release=source_params["annotation_release"],
                dir_output=self.dir_output,
            )
        elif source == "custom":
            # use already dowloaded files
            self.region_generator = CustomGenomicRegionGenerator(
                annotation_file=source_params["file_annotation"],
                sequence_file=source_params["file_sequence"],
                files_source=source_params["files_source"],
                species=source_params["species"],
                annotation_release=source_params["annotation_release"],
                genome_assembly=source_params["genome_assembly"],
                dir_output=self.dir_output,
            )
        else:
            raise ValueError(f"Source {source} not supported!")

        ##### save annotation information #####
        self.annotation_file = self.region_generator.annotation_file
        self.sequence_file = self.region_generator.sequence_file
        self.metadata["files_source"] = self.region_generator.files_source
        self.metadata["species"] = self.region_generator.species
        self.metadata["annotation_release"] = self.region_generator.annotation_release
        self.metadata["genome_assembly"] = self.region_generator.genome_assembly

        logging.info(
            f"The following annotation files are used for GTF annotation of regions: {self.annotation_file} and for fasta sequence file: {self.sequence_file} ."
        )
        logging.info(
            f"The annotations are from {self.region_generator.files_source} source, for the species: {self.region_generator.species}, release number: {self.region_generator.annotation_release} and genome assembly: {self.region_generator.genome_assembly}"
        )

    def create_oligo_database(
        self,
        regions: list,
        oligo_length_min: int,
        oligo_length_max: int,
        genomic_regions: list[
            Literal[
                "gene",
                "intergenic",
                "exon",
                "intron",
                "cds",
                "utr",
                "exon_exon_junction",
            ]
        ],
        isoform_consensus: Literal["intersection", "union"] = "union",
        min_oligos_per_region: int = 0,
        n_jobs: int = 1,
    ):
        ##### log parameters #####
        logging.info("Parameters Create Database:")
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        parameters = {i: values[i] for i in args}
        self._log_parameters(parameters)

        ##### generating the transcriptome #####
        if self.annotation_file is None or self.sequence_file is None:
            raise FileNotFoundError(
                "Annotation and Sequenec file needed to create a Transcriptome. Please use 'load_annotations()' function to provide missing files."
            )
        # length of exon_junction_size is oligo_length - 1 to continue where exons annotation ends
        fasta_files = self._parse_genomic_regions(genomic_regions=genomic_regions, block_size=oligo_length_max-1) # TODO: check the minus 1
        
        if isoform_consensus == "intersection":  # TODO: what does it mean??
            raise Exception(f"Isoform consensus: {isoform_consensus} not implemented yet.")

        ##### creating the oligo sequences #####
        oligo_sequences = OligoSequenceGenerator(dir_output=self.dir_output)
        oligo_fasta_file = oligo_sequences.create_sequences_sliding_window(
            filename_out=f"{self.log_name}_oligos",
            file_fasta_in=fasta_files,
            length_interval_sequences=(oligo_length_min, oligo_length_max),
            region_ids=regions,
        )

        ##### creating the oligo database #####
        # oligo database
        oligo_database = OligoDatabase(
            min_oligos_per_region=min_oligos_per_region,
            write_regions_with_insufficient_oligos=True,
            dir_output=self.dir_output,
        )
        # load the oligo sequences
        oligo_database.load_metadata(metadata=self.metadata)
        oligo_database.load_sequences_from_fasta(
            files_fasta=[oligo_fasta_file],
            sequence_type="oligo",
            region_ids=regions,
        )

        ##### loggig database information #####
        if self.write_removed_genes:
            logging.info(
                f"Genes with <= {min_oligos_per_region} oligos will be removed from the oligo database and their names will be stored in '{oligo_database.file_removed_regions}'."
            )

        num_genes, num_oligos = self._get_oligo_database_info(oligo_database.database)
        logging.info(
            f"Step - Generate oligos: the database contains {num_oligos} oligos from {num_genes} genes."
        )

        ##### save database #####
        if self.write_intermediate_steps:
            file_database = oligo_database.save_database(filename="oligo_database_initial.txt")
        else:
            file_database = ""

        return oligo_database, file_database

    def load_oligo_database(
        self,
        file_database: str,
        min_oligos_per_region: int = 0,
        n_jobs: int = 1,
    ):
        ##### log parameters #####
        logging.info("Parameters Load Database:")
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        parameters = {i: values[i] for i in args}
        self._log_parameters(parameters)

        ##### loading the oligo database #####
        oligo_database = OligoDatabase(
            min_oligos_per_region=min_oligos_per_region,
            dir_output=self.dir_output,
        )
        oligo_database.load_database(file_database)

        ##### loggig database information #####
        if self.write_removed_genes:
            logging.info(
                f"Genes with <= {min_oligos_per_region} oligos will be removed from the oligo database and their names will be stored in '{oligo_database.file_removed_regions}'."
            )

        num_regions, num_oligos = self._get_oligo_database_info(oligo_database.database)
        logging.info(
            f"Step - Generate oligos: the database contains {num_oligos} oligos from {num_regions} genes."
        )

        return oligo_database
    
    def _parse_genomic_regions(self, genomic_regions: dict, block_size: int = 0):
        fasta_files = []
        for genomic_region, use_flag in genomic_regions.items():
            if  use_flag: # generate the fasta file only if the flas is set to True
                if genomic_region == "gene":
                    fasta_files.append(self.region_generator.get_sequence_gene())
                elif genomic_region == "intergenic":
                    fasta_files.append(self.region_generator.get_sequence_intergenic())
                elif genomic_region == "exon":
                    fasta_files.append(self.region_generator.get_sequence_exon())
                elif genomic_region == "intron":
                    fasta_files.append(self.region_generator.get_sequence_intron())
                elif genomic_region == "cds":
                    fasta_files.append(self.region_generator.get_sequence_cds())
                elif genomic_region == "utr":
                    fasta_files.append(self.region_generator.get_sequence_utr())
                elif genomic_region == "exon_exon_junction":
                    fasta_files.append(
                        self.region_generator.get_sequence_exon_exon_junction(
                            block_size=block_size  
                        )
                    )
                else:
                    raise Exception(f"Region generator: {genomic_region} is not implemented.")
        return fasta_files
