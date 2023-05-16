import os
import sys
import yaml
import shutil
import logging
import inspect
import warnings
from pathlib import Path
from datetime import datetime

from oligo_designer_toolsuite.database import (
    EnsemblGenomicRegionGenerator,
    CustomGenomicRegionGenerator,
    NcbiGenomicRegionGenerator,
    OligoDatabase,
)


class BaseProbeDesigner:
    """This class generates all padlock probes from a transcriptome or custom file for a user-defined set of genes.
    The probe design is done in five steps:
    1. Creating all possible probes for a provided annotation and set of genes and store them in a oligo database
    2. Filter probes by a list of property filters, e.g. CG content filter
    3. Filter probes by specificity against a reference database (e.g. transcriptome)
    4. Select sets of best scoring, non-overlappign oligos for each gene
    5. Create the final ready-to-order padlock sequence

    The user can save the oligo database after each processing step and resume the pipeline later by loading an existing database into this class.
    A logger is automatically created at <output_dir>/log_padlock_probe_designer_<timestamp>.txt and logs parameters as well as number of genes/oligos after each step.

    :param dir_output: Output directory, defaults to 'output'.
    :type dir_output: str, optional
    :param write_removed_genes: write removed regions to file ``regions_with_insufficient_oligos.txt``, defaults to True
    :type write_removed_genes: bool, optional
    :param write_intermediate_steps: save oligo database after each processing step, defaults to True
    :type write_intermediate_steps: bool, optional
    """

    def __init__(
        self,
        dir_output: str = "output",
        write_removed_genes: bool = True,
        write_intermediate_steps: bool = True,
    ):
        """Constructor"""
        ##### store parameters #####
        self.dir_output = dir_output
        Path(self.dir_output).mkdir(parents=True, exist_ok=True)

        self.write_removed_genes = write_removed_genes
        self.write_intermediate_steps = write_intermediate_steps

        ##### setup logger #####
        timestamp = datetime.now()
        file_logger = os.path.join(
            self.dir_output,
            f"log_padlock_probe_designer_{timestamp.year}-{timestamp.month}-{timestamp.day}-{timestamp.hour}-{timestamp.minute}.txt",
        )
        logging.getLogger("padlock_probe_designer")
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
        self.files_source = None
        self.species = None
        self.annotation_release = None
        self.genome_assembly = None

    def _log_parameters(self, parameters):
        """Log function parameters.

        :param parameters: Dict with parameter name : parameter value pairs
        :type parameters: dict
        """
        for key, value in parameters.items():
            if key != "self":
                logging.info(f"{key} = {value}")

    def _get_probe_database_info(self, probe_database: dict):
        """Count the number of probes and genes in the database.

        :param probe_database: Database with probes.
        :type probe_database: dict
        :return: Numer of genes and probes in database.
        :rtype: int, int
        """
        genes = probe_database.keys()
        num_genes = len(genes)
        num_probes = 0
        for gene in genes:
            num_probes += len(probe_database[gene].keys())

        return num_genes, num_probes

    def _get_probe_length_min_max_from_database(self, probe_database: dict):
        """Get minimum and maximum length of probes stored in the oligo database.

        :param probe_database: Database with probes.
        :type probe_database: dict
        :return: Min and max length of probes
        :rtype: int, int
        """
        probe_length_min = sys.maxsize
        probe_length_max = 0

        for region in probe_database.keys():
            for probe in probe_database[region].keys():
                length = probe_database[region][probe]["length"]
            if length < probe_length_min:
                probe_length_min = length
            if length > probe_length_max:
                probe_length_max = length

        return probe_length_min, probe_length_max

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

        :param source: _description_
        :type source: str
        :param source_params: _description_
        :type source_params: dict
        :raises ValueError: _description_
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
        self.files_source = self.region_generator.files_source
        self.species = self.region_generator.species
        self.annotation_release = self.region_generator.annotation_release
        self.genome_assembly = self.region_generator.genome_assembly

        logging.info(
            f"The following annotation files are used for GTF annotation of regions: {self.annotation_file} and for fasta sequence file: {self.sequence_file} ."
        )
        logging.info(
            f"The annotations are from {self.files_source} source, for the species: {self.species}, release number: {self.annotation_release} and genome assembly: {self.genome_assembly}"
        )

    def create_probe_database(
        self,
        genes: list,
        probe_length_min: int,
        probe_length_max: int,
        min_probes_per_gene: int = 0,
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
        # length of exon_junction_size is probe_length - 1 to continue where exons annotation ends
        file_transcriptome = (
            self.region_generator.generate_transcript_reduced_representation(
                include_exon_junctions=True, exon_junction_size=probe_length_max
            )
        )

        ##### creating the probe database #####
        # oligo database
        probe_database = OligoDatabase(
            min_oligos_per_region=min_probes_per_gene,
            files_source=self.files_source,
            species=self.species,
            annotation_release=self.annotation_release,
            genome_assembly=self.genome_assembly,
            n_jobs=n_jobs,
            dir_output=self.dir_output,
        )
        # generate the probe sequences from gene transcripts
        probe_database.create_database(
            file_fasta=file_transcriptome,
            oligo_length_min=probe_length_min,
            oligo_length_max=probe_length_max,
            region_ids=genes,
        )

        ##### loggig database information #####
        if self.write_removed_genes:
            logging.info(
                f"Genes with <= {min_probes_per_gene} probes will be removed from the probe database and their names will be stored in '{probe_database.file_removed_regions}'."
            )

        num_genes, num_probes = self._get_probe_database_info(probe_database.database)
        logging.info(
            f"Step - Generate Probes: the database contains {num_probes} probes from {num_genes} genes."
        )

        ##### save database #####
        if self.write_intermediate_steps:
            file_database = probe_database.write_database(
                filename="probe_database_initial.txt"
            )
        else:
            file_database = ""

        return probe_database, file_database

    def load_probe_database(
        self, file_database: str, min_probes_per_gene: int = 0, n_jobs: int = 1
    ):
        ##### log parameters #####
        logging.info("Parameters Load Database:")
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        parameters = {i: values[i] for i in args}
        self._log_parameters(parameters)

        ##### loading the probe database #####
        probe_database = OligoDatabase(
            min_oligos_per_region=min_probes_per_gene,
            files_source=self.files_source,
            species=self.species,
            annotation_release=self.annotation_release,
            genome_assembly=self.genome_assembly,
            n_jobs=n_jobs,
            dir_output=self.dir_output,
        )
        probe_database.load_database(file_database)

        ##### loggig database information #####
        if self.write_removed_genes:
            logging.info(
                f"Genes with <= {min_probes_per_gene} probes will be removed from the probe database and their names will be stored in '{probe_database.file_removed_regions}'."
            )

        num_genes, num_probes = self._get_probe_database_info(probe_database.database)
        logging.info(
            f"Step - Generate Probes: the database contains {num_probes} probes from {num_genes} genes."
        )

        return probe_database
