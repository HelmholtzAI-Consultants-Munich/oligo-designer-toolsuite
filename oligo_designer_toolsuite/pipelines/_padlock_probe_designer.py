############################################
# imports
############################################


import os
import sys
import yaml
import shutil
import logging
import inspect
import warnings
from pathlib import Path
from datetime import datetime
from argparse import ArgumentParser, RawDescriptionHelpFormatter


from Bio.SeqUtils import MeltingTemp as mt

from oligo_designer_toolsuite.database import (
    CustomGenomicRegionGenerator,
    NcbiGenomicRegionGenerator,
    EnsemblGenomicRegionGenerator,
    OligoDatabase,
    ReferenceDatabase,
)
from oligo_designer_toolsuite.sequence_design import PadlockSequence
from oligo_designer_toolsuite.oligo_efficiency import (
    PadlockOligoScoring,
    PadlockSetScoring,
)
from oligo_designer_toolsuite.oligo_property_filter import (
    GCContent,
    MaskedSequences,
    MeltingTemperatureNN,
    PadlockArms,
    PropertyFilter,
)
from oligo_designer_toolsuite.oligo_selection import (
    OligosetGenerator,
    padlock_heuristic_selection,
)
from oligo_designer_toolsuite.oligo_specificity_filter import (
    Blastn,
    BowtieSeedRegion,
    ExactMatches,
    LigationRegionCreation,
    SpecificityFilter,
)
from oligo_designer_toolsuite.pipelines._padlock_probe_designer_config import (
    generate_custom_config,
    generate_ncbi_config,
    generate_ensembl_config,
)


############################################
# padlock probe design class
############################################


class PadlockProbeDesigner:
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
        logging.info("Loading Annotation Files.")
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

    def filter_probes_by_property(
        self,
        probe_database,
        GC_content_min: int = 40,
        GC_content_max: int = 60,
        Tm_min: int = 52,
        Tm_max: int = 67,
        min_arm_length: int = 10,
        max_arm_Tm_dif: int = 2,
        arm_Tm_min: int = 38,
        arm_Tm_max: int = 49,
        Tm_parameters_probe: dict = {
            "check": True,
            "strict": True,
            "c_seq": None,
            "shift": 0,
            "nn_table": "DNA_NN3",
            "tmm_table": "DNA_TMM1",
            "imm_table": "DNA_IMM1",
            "de_table": "DNA_DE1",
            "dnac1": 50,
            "dnac2": 0,
            "selfcomp": False,
            "dNTPs": 0,
            "saltcorr": 7,
            "Na": 1.25,
            "K": 75,
            "Tris": 20,
            "Mg": 10,
        },
        Tm_chem_correction_param_probe: dict = {
            "DMSO": 0,
            "DMSOfactor": 0.75,
            "fmdfactor": 0.65,
            "fmdmethod": 1,
            "GC": None,
            "fmd": 20,
        },
        n_jobs: int = 1,
    ):
        ##### log parameters #####
        logging.info("Parameters Property Filters:")
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        parameters = {i: values[i] for i in args}
        self._log_parameters(parameters)

        num_genes_before, num_probes_before = self._get_probe_database_info(
            probe_database.database
        )

        ##### preprocess melting temperature params #####
        Tm_parameters_probe["nn_table"] = getattr(mt, Tm_parameters_probe["nn_table"])
        Tm_parameters_probe["tmm_table"] = getattr(mt, Tm_parameters_probe["tmm_table"])
        Tm_parameters_probe["imm_table"] = getattr(mt, Tm_parameters_probe["imm_table"])
        Tm_parameters_probe["de_table"] = getattr(mt, Tm_parameters_probe["de_table"])

        ##### initialize the filters classes #####
        masked_sequences = MaskedSequences()
        gc_content = GCContent(
            GC_content_min=GC_content_min, GC_content_max=GC_content_max
        )
        melting_temperature = MeltingTemperatureNN(
            Tm_min=Tm_min,
            Tm_max=Tm_max,
            Tm_parameters=Tm_parameters_probe,
            Tm_chem_correction_parameters=Tm_chem_correction_param_probe,
        )
        padlock_arms = PadlockArms(
            min_arm_length=min_arm_length,
            max_arm_Tm_dif=max_arm_Tm_dif,
            arm_Tm_min=arm_Tm_min,
            arm_Tm_max=arm_Tm_max,
            Tm_parameters=Tm_parameters_probe,
            Tm_chem_correction_parameters=Tm_chem_correction_param_probe,
        )

        ##### apply property filter to the database #####
        filters = [masked_sequences, gc_content, melting_temperature, padlock_arms]
        property_filter = PropertyFilter(filters=filters)
        probe_database = property_filter.apply(
            oligo_database=probe_database, n_jobs=n_jobs
        )

        ##### save database #####
        if self.write_intermediate_steps:
            file_database = probe_database.write_database(
                filename="probe_database_property_filter.txt"
            )
        else:
            file_database = ""

        ##### loggig database information #####
        num_genes_after, num_probes_after = self._get_probe_database_info(
            probe_database.database
        )
        logging.info(
            f"Step - Filter Probes by Sequence Property: the database contains {num_probes_after} probes from {num_genes_after} genes, while {num_probes_before - num_probes_after} probes and {num_genes_before - num_genes_after} genes have been deleted in this step."
        )

        return probe_database, file_database

    def filter_probes_by_specificity(
        self,
        probe_database,
        ligation_region_size: int = 5,
        blast_word_size: int = 10,
        blast_percent_identity: int = 80,
        blast_coverage: int = 50,
        n_jobs: int = 1,
    ):
        dir_specificity = os.path.join(self.dir_output, "specificity_temporary")

        ##### log parameters #####
        logging.info("Parameters Specificity Filters:")
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        parameters = {i: values[i] for i in args}
        self._log_parameters(parameters)

        num_genes_before, num_probes_before = self._get_probe_database_info(
            probe_database.database
        )

        ##### generate transcriptome for reference #####
        # length of exon_junction_size is longer than probe length to cover bulges in alignments
        (
            probe_length_min,
            probe_length_max,
        ) = self._get_probe_length_min_max_from_database(probe_database.database)
        file_transcriptome = (
            self.region_generator.generate_transcript_reduced_representation(
                include_exon_junctions=True, exon_junction_size=2 * probe_length_max
            )
        )
        reference_database = ReferenceDatabase(
            file_fasta=file_transcriptome,
            files_source=self.region_generator.files_source,
            species=self.region_generator.species,
            annotation_release=self.region_generator.annotation_release,
            genome_assembly=self.region_generator.genome_assembly,
            dir_output=self.dir_output,
        )

        ##### intialize the filter classes #####
        exact_mathces = ExactMatches(dir_specificity=dir_specificity)
        seed_ligation = LigationRegionCreation(
            ligation_region_size=ligation_region_size
        )
        seed_region = BowtieSeedRegion(
            dir_specificity=dir_specificity,
            seed_region_creation=seed_ligation,
            strand="plus",
        )
        blastn = Blastn(
            dir_specificity=dir_specificity,
            word_size=blast_word_size,
            percent_identity=blast_percent_identity,
            coverage=blast_coverage,
            strand="plus",
        )

        ##### apply specificity filter to the database #####
        filters = [exact_mathces, blastn]  # seed_region
        specificity_filter = SpecificityFilter(filters=filters)
        probe_database = specificity_filter.apply(
            oligo_database=probe_database,
            reference_database=reference_database,
            n_jobs=n_jobs,
        )

        ##### save database #####
        if self.write_intermediate_steps:
            file_database = probe_database.write_database(
                filename="oligo_database_specificity_filters.txt"
            )
        else:
            file_database = ""

        ##### loggig database information #####
        num_genes_after, num_probes_after = self._get_probe_database_info(
            probe_database.database
        )
        logging.info(
            f"Step - Filter Probes by Specificity: the database contains {num_probes_after} probes from {num_genes_after} genes, while {num_probes_before - num_probes_after} probes and {num_genes_before - num_genes_after} genes have been deleted in this step."
        )

        shutil.rmtree(dir_specificity)

        return probe_database, file_database

    def create_probe_sets(
        self,
        probe_database,
        probeset_size_opt: int = 5,
        probeset_size_min: int = 2,
        n_sets: int = 100,
        Tm_min: int = 52,
        Tm_max: int = 67,
        Tm_opt: int = 60,
        Tm_weight: int = 1,
        GC_content_min: int = 40,
        GC_content_max: int = 50,
        GC_content_opt: int = 60,
        GC_weight: int = 1,
        n_jobs: int = 1,
    ):
        ##### log parameters #####
        logging.info("Parameters Probesets:")
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        parameters = {i: values[i] for i in args}
        self._log_parameters(parameters)

        num_genes_before, num_probes_before = self._get_probe_database_info(
            probe_database.database
        )

        ##### initialize the scoring and oligoset generator classes #####
        set_scoring = PadlockSetScoring()
        oligos_scoring = PadlockOligoScoring(
            Tm_min=Tm_min,
            Tm_opt=Tm_opt,
            Tm_max=Tm_max,
            GC_content_min=GC_content_min,
            GC_content_opt=GC_content_opt,
            GC_content_max=GC_content_max,
            Tm_weight=Tm_weight,
            GC_weight=GC_weight,
        )
        oligoset_generator = OligosetGenerator(
            oligoset_size=probeset_size_opt,
            min_oligoset_size=probeset_size_min,
            oligos_scoring=oligos_scoring,
            set_scoring=set_scoring,
            heurustic_selection=padlock_heuristic_selection,
        )

        ##### generate the oligoset #####
        probe_database = oligoset_generator.apply(
            oligo_database=probe_database, n_sets=n_sets, n_jobs=n_jobs
        )

        ##### save database #####
        if self.write_intermediate_steps:
            dir_oligosets = probe_database.write_oligosets(folder="oligosets")
            file_database = probe_database.write_database(
                filename="oligo_database_oligosets.txt"
            )
        else:
            dir_oligosets = ""
            file_database = ""

        ##### loggig database information #####
        num_genes_after, num_probes_after = self._get_probe_database_info(
            probe_database.database
        )
        logging.info(
            f"Step - Generate Oligosets: the database contains {num_probes_after} probes from {num_genes_after} genes, while {num_probes_before - num_probes_after} probes and {num_genes_before - num_genes_after} genes have been deleted in this step."
        )

        return probe_database, file_database, dir_oligosets

    def create_final_sequences(
        self,
        probe_database,
        detect_oligo_length_min: int = 18,
        detect_oligo_length_max: int = 25,
        detect_oligo_Tm_opt: int = 32,
        Tm_parameters_detection_oligo: dict = {
            "check": True,
            "strict": True,
            "c_seq": None,
            "shift": 0,
            "nn_table": "DNA_NN3",
            "tmm_table": "DNA_TMM1",
            "imm_table": "DNA_IMM1",
            "de_table": "DNA_DE1",
            "dnac1": 50,
            "dnac2": 0,
            "selfcomp": False,
            "dNTPs": 0,
            "saltcorr": 7,
            "Na": 39,
            "K": 0,
            "Tris": 0,
            "Mg": 0,
        },
        Tm_chem_correction_param_detection_oligo: dict = {
            "DMSO": 0,
            "DMSOfactor": 0.75,
            "fmdfactor": 0.65,
            "fmdmethod": 1,
            "GC": None,
            "fmd": 30,
        },
    ):
        """Generates the padlock sequences for a OligoDataset class  for which oligosets have been already computed."""

        ##### log parameters #####
        logging.info("Parameters Final Sequence Design:")
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        parameters = {i: values[i] for i in args}
        self._log_parameters(parameters)

        ##### preprocessing of the melting temperature parameters #####
        Tm_parameters_detection_oligo["nn_table"] = getattr(
            mt, Tm_parameters_detection_oligo["nn_table"]
        )
        Tm_parameters_detection_oligo["tmm_table"] = getattr(
            mt, Tm_parameters_detection_oligo["tmm_table"]
        )
        Tm_parameters_detection_oligo["imm_table"] = getattr(
            mt, Tm_parameters_detection_oligo["imm_table"]
        )
        Tm_parameters_detection_oligo["de_table"] = getattr(
            mt, Tm_parameters_detection_oligo["de_table"]
        )

        ##### initilize the padlock sequence designer class #####
        padlock_sequence = PadlockSequence(
            detect_oligo_length_min=detect_oligo_length_min,
            detect_oligo_length_max=detect_oligo_length_max,
            detect_oligo_Tm_opt=detect_oligo_Tm_opt,
            Tm_parameters=Tm_parameters_detection_oligo,
            Tm_chem_correction_parameters=Tm_chem_correction_param_detection_oligo,
            dir_output=self.dir_output,
        )

        ##### generate the final padlock sequence #####
        padlock_sequence.design_final_padlock_sequence(oligo_database=probe_database)
        logging.info(
            f"Step - Design Final Padlock Sequences: padlock sequences are stored in '{os.path.join(padlock_sequence.dir_output, 'padlock_sequences')}' directory."
        )


############################################
# commanline API
############################################


def generate_config_file(directory: str, source: str):
    directory = os.path.join(directory, "config")
    Path(directory).mkdir(parents=True, exist_ok=True)
    config_file = None
    if source == "custom":
        config_file = generate_custom_config(
            directory
        )  # function generating the config file
    elif source == "ncbi":
        config_file = generate_ncbi_config(
            directory
        )  # function generating the config file
    elif source == "ensembl":
        config_file = generate_ensembl_config(
            directory
        )  # function generating the config file
    warnings.warn(f"Default config generated automatically in '{config_file}'.")
    return config_file


def initialize_parameters(parser: ArgumentParser):
    parser.add_argument(
        "-o",
        "--output",
        help="path to the output folder, str",
        required=True,
        type=str,
        metavar="",
    )
    parser.add_argument(
        "-c",
        "--config",
        help="path to the config yaml file, str",
        default=None,
        type=str,
        metavar="",
    )
    parser.add_argument(
        "-n",
        "--n_jobs",
        help="number of cores used, int",
        type=int,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "-mog",
        "--min_oligos_per_gene",
        help="genes with less that this number of oligos are removed, int",
        type=int,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "-wg",
        "--write_removed_genes",
        help="write in a file the removed genes, bool",
        type=bool,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "-ws",
        "--write_intermediate_steps",
        help="write the oligo sequences after each step of the pipeline, bool",
        type=bool,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "-s",
        "--source",
        help="how to obtain the genomic files: download them from a server [ncbi, ensembl] or provide the files directly [custom]",
        choices=["ncbi", "ensembl", "custom"],
        default=None,
        metavar="",
    )
    parser.add_argument(
        "-tx",
        "--taxon",
        help="[source: ncbi] taxon of the species, [archaea, bacteria, fungi, invertebrate, mitochondrion, plant, plasmid, plastid, protozoa, vertebrate_mammalian, vertebrate_other, viral]",
        choices=[
            "archaea",
            "bacteria",
            "fungi",
            "invertebrate",
            "mitochondrion",
            "plant",
            "plasmid",
            "plastid",
            "protozoa",
            "vertebrate_mammalian",
            "vertebrate_other",
            "viral",
        ],
        default=None,
        metavar="",
    )
    parser.add_argument(
        "-sp",
        "--species",
        help="species name, for valid NCBI species name see https://ftp.ncbi.nlm.nih.gov/genomes/refseq/, for valid Ensembl species name see http://ftp.ensembl.org/pub/release-108/gtf/, str",
        default=None,
        metavar="",
    )
    parser.add_argument(
        "--annotation_release",
        help="release number of annotation, e.g. 'release-108' (Ensembl) or '109' (NCBI) or 'current' to use most recent annotation release, str",
        type=str,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "--genome_assembly",
        help="[source: custom] genome assembly of provided annotation, str",
        type=str,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "-fa",
        "--file_annotation",
        help="[source: custom] path to GTF file with gene annotation, str",
        type=str,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "-fs",
        "--file_sequence",
        help="[source: custom] path to FASTA file with genome sequence, str",
        type=str,
        default=None,
        metavar="",
    )

    parser.add_argument(
        "-fsrc",
        "--files_source",
        help="[source: custom] original source of the genomic files, e.g. NCBI, str",
        type=str,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "-fg",
        "--file_genes",
        help="path to file with a list of the genes that are used to generate the oligos sequences, if empty all the genes are used, str",
        type=str,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "-olm",
        "--oligo_length_min",
        help="minimum length of oligos, int",
        type=int,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "-olM",
        "--oligo_length_max",
        help="max length of oligos, int",
        type=int,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "--GC_content_min",
        help="minimum GC content of oligos, [0, 100]",
        type=float,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "--GC_content_max",
        help="maximum GC content of oligos, [0, 100]",
        type=float,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "--Tm_min",
        help="minimum melting temperature of oligos, float",
        type=float,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "--Tm_max",
        help="maximum melting temperature of oligos, float",
        type=float,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "--min_arm_length",
        help="minimum length of each arm, int",
        type=int,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "--max_arm_Tm_dif",
        help="max melting temperature difference of both arms, float",
        type=float,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "--arm_Tm_min",
        help="minimum melting temperature of each arm, float",
        type=float,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "--arm_Tm_max",
        help="max melting temperature of each arm, float",
        type=float,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "--word_size",
        help="word size for the blastn seed, int",
        type=int,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "--coverage",
        help="minimum coverage between oligos and target sequence for blastn, [0,100]",
        type=float,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "--percent_identity",
        help="maximum similarity between oligos and target sequences for blastn, [0,100]",
        type=float,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "--strand",
        help="strand of the query sequence to search, str",
        type=str,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "--ligation_region_size",
        help="size of the seed region around the ligation site for bowtie seed region filter, int",
        type=int,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "--Tm_opt",
        help="optimal melting temperature of oligos, float",
        type=float,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "--GC_content_opt",
        help="optimal GC content of oligos, [0, 100]",
        type=float,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "--Tm_weight",
        help="weight of the Tm of the oligo in the efficiency score, float",
        type=float,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "--GC_weight",
        help="weight of the GC content of the oligo in the efficiency score, float",
        type=float,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "--oligoset_size",
        help="ideal number of oligos per oligoset, int",
        type=int,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "--min_oligoset_size",
        help="minimum number of oligos per oligoset, int",
        type=int,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "--n_sets",
        help="maximum number of sets per gene, int",
        type=int,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "--detect_oligo_length_min",
        help="minimum number of oligos per oligoset, int",
        type=int,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "--detect_oligo_length_max",
        help="maximum length of detection oligo, int",
        type=int,
        default=None,
        metavar="",
    )
    parser.add_argument(
        "--detect_oligo_Tm_opt",
        help="optimal melting temperature of detection oligo, float",
        type=float,
        default=None,
        metavar="",
    )

    args = parser.parse_args()
    args = vars(args)
    if args["config"] is None:
        warnings.warn(f"No config file defined. Creating default config")
        if args["source"] is None:
            warnings.warn(f"No source was defined. Using default source: NCBI")
            args["source"] = "ncbi"
        args["config"] = generate_config_file(args["output"], args["source"])

    # read the config file
    with open(args["config"], "r") as handle:
        config = yaml.safe_load(handle)

    # update the config file values with the given one in the command line
    for param, value in args.items():
        if value is not None and param != "config":
            # overwrite the config file
            if param in [
                "file_annotation",
                "file_sequence",
                "files_source",
                "taxon",
                "species",
                "annotation_release",
                "genome_assembly",
            ]:
                # source specific parameters
                config["source_params"][param] = value
            else:
                config[param] = value

    return config


if __name__ == "__main__":
    """Command line tool to run a pipeline to design Padlock Probes, to run the tool use the command: ``padlock_probe_designer [options]``.

    The program supports two ways to recieve the input parameters:

    - command line input
    - configuration file (recommended)

    A standard configuration file can be generated using the following command ``padlock_probe_designer_config [options]``

    Since the number of input parameters requested is too high to be handled only through the command line,
    the progam will use as baseline the configuration file (it will be automatically generated if it wasn't provided)
    and the parameters specified in the command line will be overwritten with the values given.

    REMARK: melting temperature parameters can be given only through the configuration file.
    """

    # get comman line arguments
    parser = ArgumentParser(
        prog="Padlock Probe Designer",
        usage="padlock_probe_designer [options]",
        description=__doc__,
        formatter_class=RawDescriptionHelpFormatter,
    )

    config = initialize_parameters(parser)

    dir_output = os.path.abspath(config["output"])
    Path(dir_output).mkdir(parents=True, exist_ok=True)

    ##### Initialize ProbeDesigner Class #####
    probe_designer = PadlockProbeDesigner(dir_output=dir_output)

    ##### load annotations #####
    probe_designer.load_annotations(
        source=config["source"], source_params=config["source_params"]
    )

    ##### read the genes file #####
    if config["file_genes"] is None:
        warnings.warn(
            "No gene list file was provided! All genes from fasta file are used to generate the probes. This chioce can use a lot of resources."
        )
        genes = None
    else:
        with open(config["file_genes"]) as handle:
            lines = handle.readlines()
            genes = [line.rstrip() for line in lines]

    ##### create probe database #####
    probe_database, file_database = probe_designer.create_probe_database(
        genes=genes,
        probe_length_min=config["probe_length_min"],
        probe_length_max=config["probe_length_max"],
        min_probes_per_gene=config["min_probes_per_gene"],
        n_jobs=config["n_jobs"],
    )

    ##### filter probes by property #####
    probe_database, file_database = probe_designer.filter_probes_by_property(
        probe_database,
        GC_content_min=config["GC_content_min"],
        GC_content_max=config["GC_content_max"],
        Tm_min=config["Tm_min"],
        Tm_max=config["Tm_max"],
        min_arm_length=config["min_arm_length"],
        max_arm_Tm_dif=config["max_arm_Tm_dif"],
        arm_Tm_min=config["arm_Tm_min"],
        arm_Tm_max=config["arm_Tm_max"],
        Tm_parameters_probe=config["Tm_parameters_probe"],
        Tm_chem_correction_param_probe=config["Tm_chem_correction_param_probe"],
        n_jobs=config["n_jobs"],
    )

    ##### filter probes by specificity #####
    probe_database, file_database = probe_designer.filter_probes_by_specificity(
        probe_database,
        ligation_region_size=config["ligation_region_size"],
        blast_word_size=config["blast_word_size"],
        blast_percent_identity=config["blast_percent_identity"],
        blast_coverage=config["blast_coverage"],
        n_jobs=config["n_jobs"],
    )

    ##### create probe sets #####
    probe_database, file_database, dir_oligosets = probe_designer.create_probe_sets(
        probe_database,
        probeset_size_opt=config["probeset_size_opt"],
        probeset_size_min=config["probeset_size_min"],
        n_sets=config["n_sets"],
        Tm_min=config["Tm_min"],
        Tm_max=config["Tm_max"],
        Tm_opt=config["Tm_opt"],
        Tm_weight=config["Tm_weight"],
        GC_content_min=config["GC_content_min"],
        GC_content_max=config["GC_content_max"],
        GC_content_opt=config["GC_content_opt"],
        GC_weight=config["GC_weight"],
        n_jobs=config["n_jobs"],
    )

    ##### create final padlock sequence #####
    probe_designer.create_final_sequences(
        probe_database,
        detect_oligo_length_min=config["detect_oligo_length_min"],
        detect_oligo_length_max=config["detect_oligo_length_max"],
        detect_oligo_Tm_opt=config["detect_oligo_Tm_opt"],
        Tm_parameters_detection_oligo=config["Tm_parameters_detection_oligo"],
        Tm_chem_correction_param_detection_oligo=config[
            "Tm_chem_correction_param_detection_oligo"
        ],
    )
