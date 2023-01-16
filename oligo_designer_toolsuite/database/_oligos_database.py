import os
import warnings
from pathlib import Path

import pyfaidx
from joblib import cpu_count

from ..utils import FtpLoaderEnsembl, FTPLoaderNCBI, _data_parser
from ._oligos_generator import OligosGenerator
from ._transcript_generator import TranscriptGenerator


class CustomOligoDB:
    """This class generates all possible guides that can be designed for a given list of genes,
    based on the transcriptome annotation or the gene CDS or the whole genome.
    Moreover, it can write to a file and read the database containing all the oligos sequences.
    The gtf and a fasta files are passed in input.

    Species, genome_assembly, annotation_release are set to 'unknown' if thay are not given in input.

    :param oligo_length_min: minimum length of the oligos created
    :type oligo_length_min: int
    :param oligo_length_max: maximum length of the oligos created
    :type oligo_length_max: int
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
    :type n_jobs: int, optional
    :param dir_output: directory name where the results will be written
    :type dir_output: str
    :param min_oligos_per_gene: minimum number of oligos that a gene must have before it is removed from the database, defaults to 0
    :type min_oligos_per_gene: int, optional

    """

    def __init__(
        self,
        oligo_length_min: int,
        oligo_length_max: int,
        species: str = None,
        genome_assembly: str = None,
        annotation_release: str = None,
        annotation_source: str = None,
        file_annotation: str = None,
        file_sequence: str = None,
        n_jobs: int = None,
        dir_output: str = "output",
        min_oligos_per_gene: int = 0,
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
        # Initialize the file for genes with insufficient oligos
        self.file_removed_genes = os.path.join(
            self.dir_output, "genes_with_insufficient_oligos.txt"
        )
        with open(self.file_removed_genes, "a") as handle:
            handle.write(f"Gene\tPipeline step\n")
        self.species = species
        self.genome_assembly = genome_assembly
        self.annotation_release = annotation_release
        self.annotation_source = annotation_source
        self.oligo_length_max = oligo_length_max
        self.oligo_length_min = oligo_length_min
        self.min_oligos_per_gene = min_oligos_per_gene
        if n_jobs is None:
            n_jobs = cpu_count()
        self.n_jobs = n_jobs

        self.file_oligos_DB_tsv = None
        self.file_oligos_DB_gtf = None
        self.file_oligos_DB_fasta = None
        self.oligos_DB = None

        self.file_sequence = file_sequence
        self.file_annotation = file_annotation
        # create index file
        pyfaidx.Fasta(self.file_sequence)
        self.gene_transcript = None
        self.oligos = OligosGenerator(
            self.oligo_length_min,
            self.oligo_length_max,
            self.file_sequence,
            self.n_jobs,
        )
        self.oligosets = (
            {}
        )  # will be used later in the gereration of non overlpping sets

    def read_oligos_DB(
        self,
        format: str,
        file_oligos_DB_tsv: str = None,
        file_oligos_DB_gtf: str = None,
        file_oligos_DB_fasta: str = None,
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

    def write_oligos_DB(self, format: str, dir_oligos_DB: str = "oligos"):
        """
        Writes the data structure ``self.oligos_DB`` in a file.
        The fromat of the file is defined by ``format``. ``dir_oligos_DB`` is the sub-diretory of ``dir_output`` where the file will be written,

        :param format: format of file to write
        :type format: {'tsv', 'gtf'}
        :param dir_oligos_DB: path of the sub-directory where to write the file, defaults to "oligos"
        :type dir_oligos_DB: str, optional
        :return: path of the file written
        :rtype: str

        """

        if dir_oligos_DB is None:
            self.dir_oligos_DB = "oligos"
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

    def write_oligosets(self, dir_oligosets: str = "oligosets"):
        """Writes the data structure ``self.oligosets`` in a series of files, each contains the oligosets for one gene and is called "{gene}_oligosets.tsv".
        The files will be stored in a subdirectory of ``self.dir_output`` named ``dir_oligosets``.

        :param dir_oligosets: subdirectory name where the files will be stored, defaults to "oligosets"
        :type dir_oligosets: str, optional
        """

        self.dir_oligosets = os.path.join(self.dir_output, dir_oligosets)
        Path(self.dir_oligosets).mkdir(parents=True, exist_ok=True)

        for gene in self.oligosets.keys():
            file = f"{gene}_oligosets.tsv"
            path = os.path.join(self.dir_oligosets, file)
            self.oligosets[gene].to_csv(path, sep="\t", index=False)

    def __generate_gene_CDS():
        """
        Creates a fasta file containing the whole transcriptome.
        """
        # should be implemented in an external class
        raise NotImplementedError

    def create_oligos_DB(
        self,
        genes: list[str] = None,
        region: str = "gene_transcript",
    ):
        """
        Creates the DB containing all the oligo sequence extracted form the given ``region`` and belonging the the specified genes. If no genes are specified then
        will be used all the genes. The database created is not written automatically to the disk, the ``write_oligos_DB`` method hes to be called separately.

        :param genes: genes for which compute the oligos, defaults to None
        :type genes: list of str, optional
        :param region: region ofrm whihc generate the oligos, it can be 'genome', 'gene_transcript', 'gene_CDS', defaults to 'gene_transcript'
        :type region: str, optional
        :param number_batchs: oligos are computes in batches of genes, defaults to 1
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
                    self.oligo_length_max - 1, self.dir_annotation, genes
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
            self.gene_transcript = TranscriptGenerator(
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

    def remove_genes_with_insufficient_oligos(
        self, pipeline_step: str, write: bool = True
    ):
        """Deletes from the ``oligo_DB`` the genes which have less than ``min_oligos_per_gene`` oligos,
        and optionally writes them in a file with the name of the step of the pipeline at which they have been deleted.

        :param pipeline_step: name of the step of the pipeline
        :type pipeline_step: str
        """

        genes = list(self.oligos_DB.keys())
        for gene in genes:
            if len(list(self.oligos_DB[gene].keys())) <= self.min_oligos_per_gene:
                del self.oligos_DB[gene]
                if gene in self.oligosets:
                    del self.oligosets[gene]
                if write:
                    with open(self.file_removed_genes, "a") as hanlde:
                        hanlde.write(f"{gene}\t{pipeline_step}\n")


class NcbiOligoDB(CustomOligoDB):
    """
    Class to create the oligos DB using gtf and fasta files taken from the NCBI server.

    Sets species, genome_assembly, annotation_release to a predefined value if thay are not given in input.
    Dowloads the fasta and annotation files from the NCBI server and stores them in the dir_output folder.
    Moreover, it can write to a file and read the database containing all the oligos sequences.

    :param oligo_length_min: minimum length of the oligos created
    :type oligo_length_min: int
    :param oligo_length_max: maximum length of the oligos created
    :type oligo_length_max: int
    :param species: species of the files to dowload, defaults to None
    :type species: str, optional
    :param annotation_release: annotation_release of the files to dowload, defaults to None
    :type annotation_release: str, optional
    :param n_jobs: standard nr of cores used in the pipeline, if None all the available cores are used, defaults to None
    :type n_jobs: int, optional
    :param dir_output: directory where the files are saved, defaults to './output/annotation'
    :type dir_output: str, optional
    :param min_oligos_per_gene: minimum number of oligos that a gene must have before it is removed from the database, defaults to 0
    :type min_oligos_per_gene: int, optional
    """

    def __init__(
        self,
        oligo_length_min: int,
        oligo_length_max: int,
        species: str = None,
        annotation_release: str = None,
        n_jobs: int = None,
        dir_output: str = "output",
        min_oligos_per_gene: int = 0,
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

        annotation_source = "NCBI"
        dir_annotation = os.path.join(dir_output, "annotation")
        Path(dir_annotation).mkdir(parents=True, exist_ok=True)

        ftp = FTPLoaderNCBI(dir_annotation, species, annotation_release)
        file_annotation, annotation_release, genome_assembly = ftp.download_files("gtf")
        file_sequence, _, _ = ftp.download_files("fasta")

        super().__init__(
            oligo_length_min,
            oligo_length_max,
            species,
            genome_assembly,
            annotation_release,
            annotation_source,
            file_annotation,
            file_sequence,
            n_jobs,
            dir_output,
            min_oligos_per_gene,
        )


class EnsemblOligoDB(CustomOligoDB):
    """
    Class to create the oligos DB using gtf and fasta files taken from the Ensembl server.

    Sets species, genome_assembly, annotation_release to a predefined value if thay are not given in input.
    Dowloads the fasta and annotation files from the Ensemble server and stores them in the dir_output folder.
    Moreover, it can write to a file and read the database containing all the oligos sequences.

    :param oligo_length_min: minimum length of the oligos created
    :type oligo_length_min: int
    :param oligo_length_max: maximum length of the oligos created
    :type oligo_length_max: int
    :param species: species of the files to dowload, defaults to None
    :type species: str, optional
    :param genome_assembly: genome_assembly of the files to dowload, defaults to None
    :type genome_assembly: str, optional
    :param annotation_release: annotation_release of the files to dowload, defaults to None
    :type annotation_release: str, optional
    :param n_jobs: standard nr of cores used in the pipeline, if None all the available cores are used, defaults to None
    :type n_jobs: int, optional
    :param dir_output: directory where the files are saved, defaults to 'output'
    :type dir_output: str, optional
    :param min_oligos_per_gene: minimum number of oligos that a gene must have before it is removed from the database, defaults to 0
    :type min_oligos_per_gene: int, optional
    """

    def __init__(
        self,
        oligo_length_min: int,
        oligo_length_max: int,
        species: str = None,
        annotation_release: str = None,
        n_jobs: int = None,
        dir_output: str = "output",
        min_oligos_per_gene: int = 0,
    ):
        """
        Constructor
        """

        if species is None:  # change to some standard values for Ensemble
            species = "human"
            warnings.warn(f"No species defined. Using default species {species}!")

        if annotation_release is None:
            annotation_release = "current"
            warnings.warn(
                f"No annotation release defined. Using default release {annotation_release}!"
            )

        annotation_source = "Ensembl"
        dir_annotation = os.path.join(dir_output, "annotation")

        Path(dir_annotation).mkdir(parents=True, exist_ok=True)
        ftp = FtpLoaderEnsembl(dir_annotation, species, annotation_release)
        file_annotation, annotation_release, genome_assembly = ftp.download_files("gtf")
        file_sequence, _, _ = ftp.download_files("fasta")

        super().__init__(
            oligo_length_min,
            oligo_length_max,
            species,
            genome_assembly,
            annotation_release,
            annotation_source,
            file_annotation,
            file_sequence,
            n_jobs,
            dir_output,
            min_oligos_per_gene,
        )
