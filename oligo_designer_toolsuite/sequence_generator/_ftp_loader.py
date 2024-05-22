############################################
# imports
############################################

import os
import re
import gzip
import shutil
import itertools
import pandas as pd

from ftplib import FTP, error_perm
from pathlib import Path
from typing import Literal, get_args

from Bio import SeqIO

_TYPES_FILE = Literal["gff", "gtf", "fasta"]
_TYPES_SEQ = Literal["dna", "ncrna"]

############################################
# FTP Classes
############################################


class BaseFtpLoader:
    """Base class for FTP loaders.

    This class serves as a base for implementing FTP loaders, providing common functionality for
    downloading files and managing the output directory.

    :param dir_output: The directory path where files will be downloaded.
    :type dir_output: str
    """

    def __init__(self, dir_output: str):
        """Constructor for the BaseFtpLoader class."""
        self.dir_output = dir_output

    def _download(self, ftp_link: str, ftp_directory: str, file_name: str):
        """Download a file from an FTP server.

        This method connects to an FTP server, navigates to the specified directory, and downloads a file
        that matches the given filename pattern.

        :param ftp_link: The FTP link for the server.
        :type ftp_link: str
        :param ftp_directory: The directory on the FTP server where the file is located.
        :type ftp_directory: str
        :param file_name: The pattern or exact name of the file to be downloaded.
        :type file_name: str
        :return: The local path of the downloaded file, or None if the file was not found.
        :rtype: str
        """
        ftp = FTP(ftp_link)
        ftp.login()  # login to ftp server
        ftp.cwd(ftp_directory)  # move to directory

        files = ftp.nlst()
        file_output = None

        for file in files:
            if re.match(file_name, file):
                file_output = os.path.join(self.dir_output, file)
                ftp.retrbinary("RETR " + file, open(file_output, "wb").write)

        ftp.quit()

        return file_output

    def _decompress_gzip(self, file_gzip: str):
        """Decompress a gzip-compressed file.

        This method decompresses a gzip-compressed file, producing an uncompressed file in the same directory.

        :param file_gzip: The path to the gzip-compressed file.
        :type file_gzip: str
        :return: The path to the decompressed file.
        :rtype: str
        """
        file_output = file_gzip.split(".gz")[0]
        with gzip.open(file_gzip, "rb") as f_in:
            with open(file_output, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
        os.remove(file_gzip)

        return file_output

    def _download_and_decompress(self, ftp_link: str, ftp_directory: str, file_name: str):
        """Download and decompress a file from an FTP server.

        This method downloads a file from the specified FTP server, decompresses it if it is gzip-compressed,
        and returns the path to the resulting file.

        :param ftp_link: The FTP server link.
        :type ftp_link: str
        :param ftp_directory: The directory on the FTP server where the file is located.
        :type ftp_directory: str
        :param file_name: The name or pattern of the file to download.
        :type file_name: str
        :return: The path to the downloaded and decompressed file.
        :rtype: str
        """
        file_download = self._download(ftp_link, ftp_directory, file_name)
        file_unzipped = self._decompress_gzip(file_download)

        return file_unzipped

    def _check_file_type(self, file_type: _TYPES_FILE):
        """Check if the specified file type is supported.

        This method checks whether the provided file type is supported by comparing it against a predefined list
        of options.

        :param file_type: The file type to check.
        :type file_type: Literal['gff', 'gtf', 'fasta']
        :raises AssertionError: If the provided file type is not in the list of supported options.
        """
        options = get_args(_TYPES_FILE)
        assert file_type in options, f"File type not supported! '{file_type}' is not in {options}."

    def _check_sequence_nature_type(self, sequence_nature: _TYPES_SEQ):
        """Check if the provided sequence nature type is supported.

        This method checks if the provided sequence nature type is supported by comparing it to the available options.

        :param sequence_nature: The sequence nature type to be checked.
        :type sequence_nature: Literal['dna', 'ncrna']
        :raises AssertionError: If the sequence nature type is not supported.
        """
        options = get_args(_TYPES_SEQ)
        assert (
            sequence_nature in options
        ), f"Sequence nature type not supported! '{sequence_nature}' is not in {options}."


class FtpLoaderEnsembl(BaseFtpLoader):
    """FTP loader for Ensembl data.

    This class is designed to download and manage genome annotation data from Ensembl using FTP.
    It extends the functionality of the BaseFtpLoader class and provides methods for downloading
    and processing specific types of genomic files such as GFF, GTF, and FASTA.

    :param dir_output: The directory where the downloaded data will be stored.
    :type dir_output: str
    :param species: The species for which data is being downloaded (e.g., 'human', 'mouse').
    :type species: str
    :param annotation_release: The Ensembl annotation release version (e.g., '104').
    :type annotation_release: str
    """

    def __init__(self, dir_output: str, species: str, annotation_release: str):
        """Constructor for the FtpLoaderEnsembl class."""
        super().__init__(dir_output)
        self.species = species
        self.annotation_release = annotation_release
        self.assembly_name = None
        self.assembly_name_placeholder = "[^\\.]*"

        self.ftp_link = "ftp.ensembl.org"

        self.file_type_folder = {"gff": "gff3", "gtf": "gtf", "fasta": "fasta"}

        self.file_type_ending = {
            "gff": "gff3.gz",
            "gtf": "gtf.gz",
            "fasta": "dna_sm.primary_assembly.fa.gz",  # soft-masked version of the genome
        }

    def download_files(self, file_type: _TYPES_FILE, sequence_nature: _TYPES_SEQ = "dna"):
        """Download and decompress Ensembl files.

        This method downloads and decompresses Ensembl files of a specified type and sequence nature.

        :param file_type: The type of file to download.
        :type file_type: Literal['gff', 'gtf', 'fasta']
        :param sequence_nature: The sequence nature type.
        :type sequence_nature: Literal['dna', 'ncrna'], optional
        :return: Tuple containing the path to the downloaded file, annotation release, and assembly name.
        :rtype: Tuple[str, str, str]
        """
        self._check_file_type(file_type)
        self._check_sequence_nature_type(sequence_nature)

        ftp_directory, ftp_file = self._get_params(file_type, sequence_nature)
        dowloaded_file = self._download_and_decompress(self.ftp_link, ftp_directory, ftp_file)

        self.assembly_name = re.search("\\.([^\\.]*)\\.", Path(dowloaded_file).name).group().replace(".", "")

        return dowloaded_file, self.annotation_release, self.assembly_name

    def _get_params(self, file_type: _TYPES_FILE, sequence_nature: _TYPES_SEQ):
        """Get FTP parameters for downloading files.

        This method constructs the FTP directory and file name based on the provided file type and sequence nature.
        If the annotation release is set to "current," it retrieves the current Ensembl release from the README file.

        :param file_type: The type of file to be downloaded.
        :type file_type: Literal['gff', 'gtf', 'fasta']
        :param sequence_nature: The nature of the sequence.
        :type sequence_nature: Literal['dna', 'ncrna']
        :return: A tuple containing the FTP directory and file name.
        :rtype: Tuple[str, str]
        """
        Path(self.dir_output).mkdir(parents=True, exist_ok=True)

        if self.annotation_release == "current":
            file_readme = self._download(self.ftp_link, "pub/", "current_README")
            with open(file_readme, "r") as handle:
                for line in handle:
                    if line.startswith("Ensembl Release"):
                        self.annotation_release = line.strip().split(" ")[2]
            os.remove(file_readme)

        if file_type.casefold() == "fasta".casefold():
            ftp_directory = f"pub/release-{self.annotation_release}/{self.file_type_folder[file_type]}/{self.species}/{sequence_nature}/"
            if sequence_nature == "dna":
                ftp_file = f"{self.species.capitalize()}.{self.assembly_name_placeholder}.{self.file_type_ending[file_type]}"
            else:
                ftp_file = (
                    f"{self.species.capitalize()}.{self.assembly_name_placeholder}.{sequence_nature}.fa.gz"
                )
        else:
            ftp_directory = (
                f"pub/release-{self.annotation_release}/{self.file_type_folder[file_type]}/{self.species}/"
            )
            ftp_file = f"{self.species.capitalize()}.{self.assembly_name_placeholder}.{self.annotation_release}.{self.file_type_ending[file_type]}"

        return ftp_directory, ftp_file


class FtpLoaderNCBI(BaseFtpLoader):
    """A class for downloading genomic data files from the NCBI FTP server.

    The FtpLoaderNCBI class is designed to facilitate the download of genomic data files from the National Center for
    Biotechnology Information (NCBI) FTP server. It extends the functionality of the BaseFtpLoader class and provides
    methods for downloading and processing specific types of genomic files such as GFF, GTF, and FASTA.

    :param dir_output: The directory where downloaded files will be stored.
    :type dir_output: str
    :param taxon: The taxonomic identifier for the species.
    :type taxon: str
    :param species: The name of the species.
    :type species: str
    :param annotation_release: The annotation release version.
    :type annotation_release: str
    """

    def __init__(self, dir_output: str, taxon: str, species: str, annotation_release: str):
        """Constructor for the FtpLoaderNCBI class."""
        super().__init__(dir_output)
        self.taxon = taxon
        self.species = species
        self.annotation_release = annotation_release
        self.assembly_name = None
        self.assembly_accession = None

        self.ftp_link = "ftp.ncbi.nlm.nih.gov"

        self.file_type_ending = {
            "gff": "genomic.gff.gz",
            "gtf": "genomic.gtf.gz",
            "fasta": "genomic.fna.gz",  # soft-masked version of the genome
        }

        self.file_type_function = {
            "gff": self._map_chr_names_gene_annotation,
            "gtf": self._map_chr_names_gene_annotation,
            "fasta": self._map_chr_names_genome_sequence,
        }

    def download_files(self, file_type: _TYPES_FILE):
        """Download genomic data files from the NCBI FTP server.

        This method facilitates the download of specific types of genomic data files from the NCBI FTP server. It retrieves
        the necessary parameters for the specified file type, downloads the corresponding files, and performs any required
        post-processing, such as mapping chromosome names.

        :param file_type: The type of file to be downloaded.
        :type file_type: Literal['gff', 'gtf', 'fasta']
        :return: Tuple containing the path to the downloaded file, annotation release, and assembly name.
        :rtype: Tuple[str, str, str]
        """
        self._check_file_type(file_type)

        ftp_directory, ftp_file, ftp_file_chr_mapping = self._get_params(file_type)

        mapping = self._download_mapping_chr_names(ftp_directory, ftp_file_chr_mapping)
        dowloaded_file = self._download_and_decompress(self.ftp_link, ftp_directory, ftp_file)

        self.file_type_function[file_type](dowloaded_file, mapping)

        return dowloaded_file, self.annotation_release, self.assembly_name

    def _get_params(self, file_type: _TYPES_FILE):
        """Get FTP parameters for downloading genomic data files.

        This method retrieves the FTP parameters necessary for downloading specific types of genomic data files from the
        NCBI FTP server. It checks the file type, creates the necessary directories, and determines the appropriate FTP
        directory, file names, and file paths.

        :param file_type: The type of file to be downloaded.
        :type file_type: Literal['gff', 'gtf', 'fasta']
        :return: A tuple containing the FTP directory, file name for the genomic data file, and file name for the chromosome name mapping file.
        :rtype: Tuple[str, str, str]
        """
        Path(self.dir_output).mkdir(parents=True, exist_ok=True)

        ftp_directory = "genomes/refseq/" + self.taxon + "/" + self.species + "/annotation_releases/"

        if self.annotation_release == "current":
            ftp_directory = ftp_directory + "current/"
            ftp = FTP(
                self.ftp_link
            )  # inside current dir there is the directory containing the annotation release
            ftp.login()  # login to ftp server
            ftp.cwd(ftp_directory)  # move to directory
            self.annotation_release = ftp.nlst()[0]
            ftp.quit()

        ftp_directory = ftp_directory + f"{self.annotation_release}/"

        file_readme = self._download(self.ftp_link, ftp_directory, f"README_.*{self.annotation_release}")
        with open(file_readme, "r") as handle:
            for line in handle:
                if line.startswith("ASSEMBLY NAME:"):
                    self.assembly_name = line.strip().split("\t")[1]
                if line.startswith("ASSEMBLY ACCESSION:"):
                    self.assembly_accession = line.strip().split("\t")[1]
                    break
        os.remove(file_readme)

        # we need this check here because NCBI changed the folder structure for releases > 110
        ftp = FTP(self.ftp_link)
        ftp.login()
        try:
            ftp.cwd(ftp_directory + f"{self.assembly_accession}_{self.assembly_name}")  # move to directory
            ftp_directory = ftp_directory + f"{self.assembly_accession}_{self.assembly_name}"
        except error_perm:
            ftp_directory = ftp_directory
        ftp.quit()

        ftp_file = f"{self.assembly_accession}_{self.assembly_name}_{self.file_type_ending[file_type]}"
        ftp_file_chr_mapping = f"{self.assembly_accession}_{self.assembly_name}_assembly_report.txt"

        return ftp_directory, ftp_file, ftp_file_chr_mapping

    def _download_mapping_chr_names(self, ftp_directory: str, ftp_file_chr_mapping: str):
        """Download and parse the chromosome name mapping file.

        This method downloads the chromosome name mapping file from the specified FTP directory, parses the file to extract
        the relevant information, and returns a dictionary mapping RefSeq accessions to chromosome names.

        :param ftp_directory: The FTP directory containing the chromosome name mapping file.
        :type ftp_directory: str
        :param ftp_file_chr_mapping: The name of the chromosome name mapping file.
        :type ftp_file_chr_mapping: str
        :return: A dictionary mapping RefSeq accessions to chromosome names.
        :rtype: Dict
        """
        file_mapping = self._download(self.ftp_link, ftp_directory, ftp_file_chr_mapping)

        # skip comment lines but keep last comment line for header
        with open(file_mapping) as handle:
            *_comments, names = itertools.takewhile(lambda line: line.startswith("#"), handle)
            names = names[1:].split()

        assembly_report = pd.read_table(file_mapping, names=names, sep="\t", comment="#")

        mapping_chromosome = assembly_report[assembly_report["Sequence-Role"] == "assembled-molecule"]
        mapping_chromosome = pd.Series(
            mapping_chromosome["Sequence-Name"].values,
            index=mapping_chromosome["RefSeq-Accn"],
        ).to_dict()

        mapping_scaffolds = assembly_report[assembly_report["Sequence-Role"] != "assembled-molecule"]
        mapping_scaffolds = pd.Series(
            mapping_scaffolds["GenBank-Accn"].values,
            index=mapping_scaffolds["RefSeq-Accn"],
        ).to_dict()

        mapping = mapping_chromosome
        mapping.update(mapping_scaffolds)

        return mapping

    def _map_chr_names_gene_annotation(self, ftp_file: str, mapping: dict):
        """Map chromosome names in a gene annotation file.

        This method reads a gene annotation file in GTF format, maps chromosome names using the provided mapping, and
        writes the modified annotation to the same file.

        :param ftp_file: The path to the gene annotation file.
        :type ftp_file: str
        :param mapping: A dictionary mapping RefSeq accessions to chromosome names.
        :type mapping: Dict
        """
        file_tmp = os.path.join(self.dir_output, "temp.gtf")

        # write comment lines to new file
        with open(file_tmp, "w") as handle_out:
            with open(ftp_file) as handle_in:
                *_comments, names = itertools.takewhile(lambda line: line.startswith("#"), handle_in)
                handle_out.write(names)

            # read gtf file without comment lines
            gene_annotation = pd.read_table(
                ftp_file,
                names=[
                    "seqid",
                    "source",
                    "type",
                    "start",
                    "end",
                    "score",
                    "strand",
                    "phase",
                    "attributes",
                ],
                sep="\t",
                comment="#",
            )

            # replace ncbi with genbank chromosome annotation
            gene_annotation["seqid"] = gene_annotation["seqid"].map(mapping)
            gene_annotation.dropna(inplace=True)  # drop if no mapping exists

            gene_annotation.to_csv(handle_out, sep="\t", header=False, index=False)
        os.replace(file_tmp, ftp_file)

    def _map_chr_names_genome_sequence(self, ftp_file: str, mapping: dict):
        """Map chromosome names in a genome sequence file.

        This method reads a genome sequence file in FASTA format, maps chromosome names using the provided mapping, and
        writes the modified sequence to the same file.

        :param ftp_file: The path to the genome sequence file.
        :type ftp_file: str
        :param mapping: A dictionary mapping RefSeq accessions to chromosome names.
        :type mapping: Dict
        """
        file_tmp = os.path.join(self.dir_output, "temp.fna")

        with open(file_tmp, "w") as handle:
            for chromosome_sequnece in SeqIO.parse(ftp_file, "fasta"):
                accession_number = chromosome_sequnece.id
                if accession_number in mapping:
                    chromosome_sequnece.id = mapping[accession_number]
                    chromosome_sequnece.name = mapping[accession_number]
                    chromosome_sequnece.description = chromosome_sequnece.description.replace(
                        accession_number, mapping[accession_number]
                    )
                    SeqIO.write(chromosome_sequnece, handle, "fasta")
                else:
                    self.logging.info("No mapping for accession number: {}".format(accession_number))

        os.replace(file_tmp, ftp_file)
