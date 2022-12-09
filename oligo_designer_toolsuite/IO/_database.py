import os
import shutil
import warnings
from pathlib import Path

import pyfaidx
from joblib import cpu_count

from ..oligo_transcript_generation import GeneTranscript, Oligos
from ..utils import FtpLoaderEnsembl, FTPLoaderNCBI, _data_parser


class CustomDB:
    """This class generates all possible guides that can be designed for a given list of genes,
    based on the transcriptome annotation or the gene CDS or the whole genome and the reference fasta file.
    The gtf and a fasta files are passed in input

    Sets species, genome_assembly, annotation_release to 'unknown' if thay are not given in input. Saves the path of the user defined annotation and fasta file and initializes the class argumenets.

    :param probe_length_min: minimum length of the probes created
    :type probe_length_min: int
    :param probe_length_max: maximum length of the probes created
    :type probe_length_max: int
    :param species: species of the fasta and gtf files, defaults to None
    :type species: str, optional
    :param genome_assembly: genome_assembly of the fasta and gtf files, defaults to None
    :type genome_assembly: str, optional
    :param annotation_release: annotation_release of the fasta and gtf files, defaults to None
    :type annotation_release: str, optional
    :param annotation_source: source of the fasta and gtf files, defaults to None
    :type annotation_source: str, optional
    :param file_annotation: path to the gtf annotation file, defaults to None
    :type file_annotation: str, optional
    :param file_sequence: path to the fasta file, defaults to None
    :type file_sequence: str, optional
    :param n_jobs: standard nr of cores used in the pipeline, if None all the available cores are used, defaults to None
    :type n_jobs: int
    :param dir_output: directory name where the results will be written
    :type dir_output: str

    """

    def __init__(
        self,
        probe_length_min,
        probe_length_max,
        species=None,
        genome_assembly=None,
        annotation_release=None,
        annotation_source=None,
        file_annotation=None,
        file_sequence=None,
        n_jobs=None,
        dir_output="output",
        min_probes_per_gene=0,
    ):
        """
        Constructor
        """
        if species is None:
            species = "unknown"
            warnings.warn("Species not specified.")

        if genome_assembly is None:
            genome_assembly = "unknown"
            warnings.warn(f"Genome assembly not specified.")

        if annotation_release is None:
            annotation_release = "unknown"
            warnings.warn(f"Annotation release not specified.")

        if annotation_source is None:
            annotation_source = "Custom"
            warnings.warn(f"Annotation source not specified.")

        # check the files format
        if file_annotation == None:
            raise ValueError("Annotation File not defined!")

        if file_sequence == None:
            raise ValueError("Sequence File not defined!")

        if not _data_parser.check_gtf_format(file_annotation):
            raise ValueError("Annotation File has incorrect format!")

        if not _data_parser.check_fasta_format(file_sequence):
            raise ValueError("Sequence File has incorrect format!")

        self.dir_output = dir_output
        self.dir_annotation = os.path.join(dir_output, "annotation")
        Path(self.dir_annotation).mkdir(parents=True, exist_ok=True)
        # Initialize the file for genes with insufficient probes
        self.file_removed_genes = os.path.join(
            self.dir_output, "genes_with_insufficient_probes.txt"
        )
        with open(self.file_removed_genes, "a") as handle:
            handle.write(f"Gene\tPipeline step\n")
        self.species = species
        self.genome_assembly = genome_assembly
        self.annotation_release = annotation_release
        self.annotation_source = annotation_source
        self.probe_length_max = probe_length_max
        self.probe_length_min = probe_length_min
        self.min_probes_per_gene = min_probes_per_gene
        if n_jobs is None:
            n_jobs = cpu_count()
        self.n_jobs = n_jobs

        self.file_reference_DB = None
        self.file_oligos_DB_tsv = None
        self.file_oligos_DB_gtf = None
        self.file_oligos_DB_fasta = None
        self.oligos_DB = None

        self.file_sequence = file_sequence
        self.file_annotation = file_annotation
        # create index file
        pyfaidx.Fasta(self.file_sequence)
        self.gene_transcript = None
        self.oligos = Oligos(
            self.probe_length_min,
            self.probe_length_max,
            self.file_sequence,
            self.n_jobs,
        )
        self.probesets = (
            {}
        )  # will be used later in the gereration of non overlpping sets

    def read_reference_DB(self, file_reference_DB):
        """Saves the path of a previously generated reference DB in the ``self.file_reference_DB`` attribute.

        :param file_reference_DB: path of the reference_DB file
        :type file_reference_DB: str

        """
        if os.path.exists(file_reference_DB):
            if _data_parser.check_fasta_format(file_reference_DB):
                self.file_reference_DB = file_reference_DB
            else:
                raise ValueError("Database has incorrect format!")
        else:
            raise ValueError("Database file does not exist!")

    def read_oligos_DB(
        self,
        format,
        file_oligos_DB_tsv=None,
        file_oligos_DB_gtf=None,
        file_oligos_DB_fasta=None,
    ):
        """
        Create the oligo db dictionary from a file. It can take both a tsv file of a gtf and fasta file and the format of file to process is defined by ``format``.

        :param format: format of file to process
        :type format: {'tsv', 'gtf'}
        :param file_oligos_DB_gtf: Path to the file.
        :type file_oligos_DB_gtf: str
        :param file_oligos_DB_fasta: Path to the file.
        :type file_oligos_DB_fasta: str
        :param file_oligos_DB_tsv: path of the oligos_DB file
        :type file_oligos_DB_tsv: str

        """

        if format == "tsv":
            self.oligos_DB = _data_parser.read_oligos_DB_tsv(file_oligos_DB_tsv)
            self.file_oligos_DB_tsv = (
                file_oligos_DB_tsv  # already checked if it is a tsv file
            )
        elif format == "gtf":
            self.oligos_DB = _data_parser.read_oligos_DB_gtf(
                file_oligos_DB_gtf, file_oligos_DB_fasta
            )
            self.file_oligos_DB_gtf = file_oligos_DB_gtf
            self.file_oligos_DB_fasta = file_oligos_DB_fasta
        else:
            raise ValueError(f"{format} not recognized as a format!")

    def write_oligos_DB(self, format, dir_oligos_DB=None):
        """
        Writes the data structure self.oligos_DB in a file.
        The fromat of the file is defined by ``format``. ``dir_oligos_DB`` is the sub-diretory of ``dir_output`` where the file will be written,
        if None it will be set as the ``dir_annotation``.

        :param format: format of file to write
        :type format: {'tsv', 'gtf'}
        :param dir_oligos_DB: path of the sub-directory where to write the file, defaults to None
        :type dir_oligos_DB: str, optional
        :return: path of the file written
        :rtype: str

        """

        if dir_oligos_DB is None:
            self.dir_oligos_DB = self.dir_annotation
        else:
            self.dir_oligos_DB = os.path.join(self.dir_output, dir_oligos_DB)
            Path(self.dir_oligos_DB).mkdir(parents=True, exist_ok=True)

        if format == "tsv":
            self.file_oligos_DB_tsv = os.path.join(
                self.dir_oligos_DB,
                self.file_name_oligos_DB_tsv,
            )
            _data_parser.write_oligos_DB_tsv(self.oligos_DB, self.file_oligos_DB_tsv)
            return self.file_name_oligos_DB_tsv
        elif format == "gtf":
            self.file_oligos_DB_gtf = os.path.join(
                self.dir_oligos_DB,
                self.file_name_oligos_DB_gtf,
            )
            self.file_oligos_DB_fasta = os.path.join(
                self.dir_oligos_DB,
                self.file_name_oligos_DB_fasta,
            )
            _data_parser.write_oligos_DB_gtf(
                self.oligos_DB, self.file_oligos_DB_gtf, self.file_oligos_DB_fasta
            )
            return self.file_name_oligos_DB_gtf, self.file_name_oligos_DB_fasta
        else:
            raise ValueError(f"{format} not recognized as a format!")

    def create_reference_DB(
        self,
        region="gene_transcript",
        block_size=None,
        dir_reference_DB=None,
    ):
        """
        Creates a fasta file for each of the region selected (genome, gene_transcript, gene_CDS) which will be used for alignements, default is "gene_transcript".
        If not specified the exon junctions size is set to ``probe_length_max`` + 5. ``dir_reference_DB`` is the subdirectory of dir_out where the reference file will be written,
        if None it will be set to dir_annotation.

        :param region: the region to use for the reference DB. Possible values are "genome", "gene_transcript", "gene_CDS"
        :type region: str
        :param block_size: size of the exon junctions. When specified as None, the block size is set to ``probe_length_max`` + 5, defaults to None
        :type block_size: int, optional
        :param dir_reference_DB: path of the sub-directory where to write the file, defaults to None
        :type dir_reference_DB: str, optional
        :return: path of the file written
        :rtype: str

        """
        self.gene_transcript = GeneTranscript(self.file_sequence, self.file_annotation)
        if dir_reference_DB is None:
            dir_reference_DB = self.dir_annotation
        else:
            dir_reference_DB = os.path.join(self.dir_output, dir_reference_DB)
            Path(dir_reference_DB).mkdir(parents=True, exist_ok=True)

        def get_files_fasta(region, dir_reference_DB, file_reference_DB):
            """
            Generates the fasta files that will compose the reference_DB and writes it in ``file_reference_DB``

            :param region: the region to use for the reference DB
            :type region: {'genome', 'gene_transcript', 'gene_CDS'}
            :param dir_reference_DB: path of the directory where to write the intermediate files.
            :type dir_reference_DB: str
            :param file_reference_DB: path of the file where to write the reference_DB.
            :type file_reference_DB: str

            """
            if region == "genome":
                shutil.copyfile(self.file_sequence, file_reference_DB)
            elif region == "gene_transcript":
                self.gene_transcript.generate_for_reference(
                    block_size, file_reference_DB, dir_reference_DB
                )  # call the outer class to generate the gene trascript
            elif region == "gene_CDS":
                # generate gene cds
                warnings.warn("Gene CDS not implemented yet")
            else:
                raise ValueError(
                    f"The given region does not exists. You selected {region} but only 'genome', 'gene_transcript', 'gene_CDS' are available."
                )

        self.file_reference_DB = os.path.join(
            dir_reference_DB,
            f"reference_DB_{self.species}_{self.genome_assembly}_{self.annotation_source}_release_{self.annotation_release}_{region}.fna",
        )

        if block_size is None:  # when not specified define it automatically
            block_size = self.probe_length_max + 5

        get_files_fasta(region, dir_reference_DB, self.file_reference_DB)
        return self.file_reference_DB

    def __generate_gene_CDS():
        """
        Creates a fasta file containing the whole transcriptome.
        """
        # should be implemented in an external class
        raise NotImplementedError

    def create_oligos_DB(
        self,
        genes=None,
        region="gene_transcript",
    ):
        """
        Creates the DB containing all the oligo sequence extracted form the given ``region`` and belonging the the specified genes. If no genes are specified then
        will be used all the genes. The database created is not written automatically to the disk, the ``write_oligos_DB`` method hes to be called separately.

        :param genes: genes for which compute the probes, defaults to None
        :type genes: list of str, optional
        :param region: region ofrm whihc generate the probes, it can be 'genome', 'gene_transcript', 'gene_CDS', defaults to 'gene_transcript'
        :type region: str, optional
        :param number_batchs: probes are computes in batches of genes, defaults to 1
        :type number_batchs: int, optional

        """

        def create_target_region(region, genes):
            """cretares the annotation file for the specified region."""

            if region == "genome":
                file_region_annotation = None
                warnings.warn("Genome not implemented yet")
                # TODO
            elif region == "gene_transcript":
                file_region_annotation = self.gene_transcript.generate_for_oligos(
                    self.probe_length_max - 1, self.dir_annotation, genes
                )
            elif region == "gene_CDS":
                file_region_annotation = None
                warnings.warn("Gene CDS not implemented yet")
                # TODO
            else:
                raise ValueError(
                    f"The given region does not exists. You selected {region}"
                )
            return file_region_annotation

        self.file_name_oligos_DB_tsv = f"oligo_DB_{self.species}_{self.genome_assembly}_{self.annotation_source}_release_{self.annotation_release}_{region}.tsv"
        self.file_name_oligos_DB_gtf = f"oligo_DB_{self.species}_{self.genome_assembly}_{self.annotation_source}_release_{self.annotation_release}_{region}.gtf"
        self.file_name_oligos_DB_fasta = f"oligo_DB_{self.species}_{self.genome_assembly}_{self.annotation_source}_release_{self.annotation_release}_{region}.fasta"

        if self.gene_transcript is None:
            self.gene_transcript = GeneTranscript(
                self.file_sequence, self.file_annotation
            )
        file_region_annotation = create_target_region(region, genes)

        if genes is None:
            genes = self.gene_transcript.get_genes_from_annotation()

        self.oligos_DB = self.oligos.generate(
            file_region_annotation, genes, self.dir_annotation
        )
        # clean folder
        os.remove(file_region_annotation)

    def remove_genes_with_insufficient_probes(self, pipeline_step, write=True):
        """Deletes from the ``oligo_DB`` the genes which have less than ``min_probes_per_gene`` probes,
        and optionally writes them in a file with the name of the step of the pipeline at which they have been deleted.

        :param pipeline_step: name of the step of the pipeline
        :type pipeline_step: str
        """

        genes = list(self.oligos_DB.keys())
        for gene in genes:
            if len(list(self.oligos_DB[gene].keys())) <= self.min_probes_per_gene:
                del self.oligos_DB[gene]
                if write:
                    with open(self.file_removed_genes, "a") as hanlde:
                        hanlde.write(f"{gene}\t{pipeline_step}\n")


class NcbiDB(CustomDB):
    """
    Class to create reference and oligos DB using gtf and fasta files taken from the NCBI server.

    Sets species, genome_assembly, annotation_release to a predefined value if thay are not given in input.
    Dowloads the fasta and annotation files from the NCBI server and stores them in the dir_output folder

    :param probe_length_min: minimum length of the probes created
    :type probe_length_min: int
    :param probe_length_max: maximum length of the probes created
    :type probe_length_max: int
    :param species: species of the files to dowload, defaults to None
    :type species: str, optional
    :param annotation_release: annotation_release of the files to dowload, defaults to None
    :type annotation_release: str, optional
    :param dir_output: directory where the files are saved, defaults to './output/annotation'
    :type dir_output: str, optional
    :param filters: list of filters classes already initialized, defaults to None
    :type filters: list of classes, optional
    """

    def __init__(
        self,
        probe_length_min,
        probe_length_max,
        filters,
        species=None,
        annotation_release=None,
        dir_output="output",
    ):
        """
        Constructor
        """
        if species is None:
            species = "human"
            warnings.warn(f"No species defined. Using default species {species}!")

        if annotation_release is None:
            annotation_release = "current"
            warnings.warn(
                f"No annotation release defined. Using default release {annotation_release}!"
            )

        genome_assembly = "GRCh38"
        annotation_source = "NCBI"
        dir_annotation = os.path.join(dir_output, "annotation")
        Path(dir_annotation).mkdir(parents=True, exist_ok=True)

        ftp = FTPLoaderNCBI(dir_annotation, species, annotation_release)
        file_annotation = ftp.download_files("gtf")
        file_sequence = ftp.download_files("fasta")

        super().__init__(
            probe_length_min,
            probe_length_max,
            filters,
            species,
            genome_assembly,
            annotation_release,
            annotation_source,
            file_annotation,
            file_sequence,
            dir_output,
        )


class EnsemblDB(CustomDB):
    """
    Class to create reference and oligos DB using gtf and fasta files taken from the Ensembl server.

    Sets species, genome_assembly, annotation_release to a predefined value if thay are not given in input.
    Dowloads the fasta and annotation files from the Ensemble server and stores them in the dir_output folder.

    :param probe_length_min: minimum length of the probes created
    :type probe_length_min: int
    :param probe_length_max: maximum length of the probes created
    :type probe_length_max: int
    :param species: species of the files to dowload, defaults to None
    :type species: str, optional
    :param genome_assembly: genome_assembly of the files to dowload, defaults to None
    :type genome_assembly: str, optional
    :param annotation_release: annotation_release of the files to dowload, defaults to None
    :type annotation_release: str, optional
    :param dir_output: directory where the files are saved, defaults to 'output'
    :type dir_output: str, optional
    :param filters: list of filters classes already initialized, defaults to None
    :type filters: list of classes, optional
    """

    def __init__(
        self,
        probe_length_min,
        probe_length_max,
        filters,
        species=None,
        genome_assembly=None,
        annotation_release=None,
        dir_output="output",
    ):
        """
        Constructor
        """

        if species is None:  # change to some standard values for Ensemble
            species = "human"
            warnings.warn(f"No species defined. Using default species {species}!")

        if genome_assembly is None:
            genome_assembly = "GRCh38"
            warnings.warn(
                f"No genome assembly defined. Using default assembly {genome_assembly}!"
            )

        if annotation_release is None:
            annotation_release = "current"
            warnings.warn(
                f"No annotation release defined. Using default release {annotation_release}!"
            )

        annotation_source = "Ensembl"
        dir_annotation = os.path.join(dir_output, "annotation")

        Path(dir_annotation).mkdir(parents=True, exist_ok=True)
        ftp = FtpLoaderEnsembl(
            dir_annotation, species, genome_assembly, annotation_release
        )
        file_annotation = ftp.download_files("gtf")
        file_sequence = ftp.download_files("fasta")

        super().__init__(
            probe_length_min,
            probe_length_max,
            filters,
            species,
            genome_assembly,
            annotation_release,
            annotation_source,
            file_annotation,
            file_sequence,
            dir_output,
        )
