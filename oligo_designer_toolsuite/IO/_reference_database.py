import os
import shutil
import warnings
from pathlib import Path

import pyfaidx

from ..oligo_transcript_generation import GeneTranscript
from ..utils import FtpLoaderEnsembl, FTPLoaderNCBI, _data_parser


class CustomReferenceDB:
    """This class generates reference fasta file used for the alignement methods.
    The gtf and a fasta files are passed in input.

    Sets species, genome_assembly, annotation_release to 'unknown' if thay are not given in input.

    :param file_annotation: path to the gtf annotation file, defaults to None
    :type file_annotation: str, optional
    :param file_sequence: path to the fasta file, defaults to None
    :type file_sequence: str, optional
    :param species: species of the fasta and gtf files, defaults to None
    :type species: str, optional
    :param genome_assembly: genome_assembly of the fasta and gtf files, defaults to None
    :type genome_assembly: str, optional
    :param annotation_release: annotation_release of the fasta and gtf files, defaults to None
    :type annotation_release: str, optional
    :param annotation_source: source of the fasta and gtf files, defaults to None
    :type annotation_source: str, optional
    :param dir_output: directory name where the results will be written
    :type dir_output: str

    """

    def __init__(
        self,
        file_annotation: str = None,
        file_sequence: str = None,
        species: str = None,
        genome_assembly: str = None,
        annotation_release: str = None,
        annotation_source: str = None,
        dir_output: str = "output",
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
        self.species = species
        self.genome_assembly = genome_assembly
        self.annotation_release = annotation_release
        self.annotation_source = annotation_source

        self.file_reference_DB = None

        self.file_sequence = file_sequence
        self.file_annotation = file_annotation
        # create index file
        pyfaidx.Fasta(self.file_sequence)
        self.gene_transcript = None

    def read_reference_DB(self, file_reference_DB: str):
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

    def create_reference_DB(
        self,
        block_size: int,
        region: str = "gene_transcript",
        dir_reference_DB: str = "reference",
    ):
        """
        Creates a fasta file for each of the region selected (genome, gene_transcript, gene_CDS) which will be used for alignements, default is "gene_transcript".
        ``dir_reference_DB`` is the subdirectory of dir_out where the reference file will be written.
        For the "gene_transcript", the block size are the additional bases added on the estremes of each exon.

        :param region: the region to use for the reference DB. Possible values are "genome", "gene_transcript", "gene_CDS"
        :type region: str
        :param block_size: size of the exon junctions.
        :type block_size: int
        :param dir_reference_DB: path of the sub-directory where to write the file, defaults to "reference"
        :type dir_reference_DB: str, optional
        :return: path of the file written
        :rtype: str

        """
        if self.gene_transcript is None:
            self.gene_transcript = GeneTranscript(
                self.file_sequence, self.file_annotation
            )

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

        get_files_fasta(region, dir_reference_DB, self.file_reference_DB)
        return self.file_reference_DB

    def __generate_gene_CDS():
        """
        Creates a fasta file containing the whole transcriptome.
        """
        # should be implemented in an external class
        raise NotImplementedError


class NcbiReferenceDB(CustomReferenceDB):
    """
    Class to create reference DB using gtf and fasta files taken from the NCBI server.

    Sets species, genome_assembly, annotation_release to a predefined value if thay are not given in input.
    Dowloads the fasta and annotation files from the NCBI server and stores them in the dir_output folder

    :param species: species of the files to dowload, defaults to None
    :type species: str, optional
    :param annotation_release: annotation_release of the files to dowload, defaults to None
    :type annotation_release: str, optional
    :param dir_output: directory where the files are saved, defaults to './output/annotation'
    :type dir_output: str, optional
    """

    def __init__(
        self,
        species: str = None,
        annotation_release: str = None,
        dir_output: str = "output",
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
            file_annotation,
            file_sequence,
            species,
            genome_assembly,
            annotation_release,
            annotation_source,
            dir_output,
        )


class EnsemblReferenceDB(CustomReferenceDB):
    """
    Class to create reference DB using gtf and fasta files taken from the Ensembl server.

    Sets species, genome_assembly, annotation_release to a predefined value if thay are not given in input.
    Dowloads the fasta and annotation files from the Ensemble server and stores them in the dir_output folder.

    :param species: species of the files to dowload, defaults to None
    :type species: str, optional
    :param annotation_release: annotation_release of the files to dowload, defaults to None
    :type annotation_release: str, optional
    :param dir_output: directory where the files are saved, defaults to './output/annotation'
    :type dir_output: str, optional
    """

    def __init__(
        self,
        species: str = None,
        annotation_release: str = None,
        dir_output: str = "output",
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
            file_annotation,
            file_sequence,
            species,
            genome_assembly,
            annotation_release,
            annotation_source,
            dir_output,
        )
