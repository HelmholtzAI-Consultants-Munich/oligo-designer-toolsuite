############################################
# imports
############################################

import gzip
import itertools
import os
import re
import shutil
from ftplib import FTP, error_perm
from pathlib import Path
from typing import Tuple, get_args

import pandas as pd
from Bio import SeqIO

from oligo_designer_toolsuite._constants import _TYPES_FILE, _TYPES_FILE_SEQ

############################################
# FTP Classes
############################################


class BaseFtpLoader:
    """
    A base class for downloading files via FTP and postprocessing the downloaded files.

    :param dir_output: The directory path where the downloaded files will be saved.
    :type dir_output: str
    """

    def __init__(self, dir_output: str) -> None:
        """Constructor for the BaseFtpLoader class."""
        self.dir_output = dir_output

    def _download(self, ftp_link: str, ftp_directory: str, file_name: str) -> str:
        """
        Downloads a file from an FTP server.

        :param ftp_link: The link to the FTP server.
        :type ftp_link: str
        :param ftp_directory: The directory on the FTP server where the file is located.
        :type ftp_directory: str
        :param file_name: The name of the file to download, can be a regex pattern.
        :type file_name: str
        :return: The path to the downloaded file.
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

    def _decompress_gzip(self, file_gzip: str) -> str:
        """
        Decompresses a .gz file and removes the original compressed file.

        :param file_gzip: The path to the .gz file to be decompressed.
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

    def _download_and_decompress(self, ftp_link: str, ftp_directory: str, file_name: str) -> str:
        """
        Downloads a file from an FTP server, decompresses it, and returns the path to the decompressed file.

        :param ftp_link: The link to the FTP server.
        :type ftp_link: str
        :param ftp_directory: The directory on the FTP server where the file is located.
        :type ftp_directory: str
        :param file_name: The name of the file to download, can be a regex pattern.
        :type file_name: str
        :return: The path to the decompressed file.
        :rtype: str
        """
        file_download = self._download(ftp_link, ftp_directory, file_name)
        file_unzipped = self._decompress_gzip(file_download)

        return file_unzipped

    def _check_file_type(self, file_type: _TYPES_FILE) -> None:
        """
        Checks if the provided file type is supported.

        :param file_type: The type of file to check.
        :type file_type: _TYPES_FILE ["gff", "gtf", "fasta"]
        """
        options = get_args(_TYPES_FILE)
        assert file_type in options, f"File type not supported! '{file_type}' is not in {options}."

    def _check_sequence_nature_type(self, sequence_nature: _TYPES_FILE_SEQ) -> None:
        """
        Checks if the provided sequence nature type is supported.

        :param sequence_nature: The type of sequence nature to check.
        :type sequence_nature: _TYPES_FILE_SEQ["dna", "ncrna"]
        """
        options = get_args(_TYPES_FILE_SEQ)
        assert (
            sequence_nature in options
        ), f"Sequence nature type not supported! '{sequence_nature}' is not in {options}."


class FtpLoaderEnsembl(BaseFtpLoader):
    """
    A class for downloading genomic data from the Ensembl FTP server.

    The `FtpLoaderEnsembl` class is designed to facilitate the retrieval of genomic data, such as GFF, GTF, and FASTA files,
    for a specific species and annotation release from the Ensembl FTP server. The class handles the construction of FTP paths and file names,
    and manages the download and decompression of files.

    :param dir_output: The directory where the downloaded files will be saved.
    :type dir_output: str
    :param species: The species for which the genomic data is to be downloaded (e.g., 'human', 'mouse').
    :type species: str
    :param annotation_release: The Ensembl annotation release version (e.g., '104').
    :type annotation_release: str
    """

    def __init__(self, dir_output: str, species: str, annotation_release: str) -> None:
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

    def download_files(
        self, file_type: _TYPES_FILE, sequence_nature: _TYPES_FILE_SEQ = "dna"
    ) -> Tuple[str, str, str]:
        """
        Downloads and decompresses genomic data files from the Ensembl FTP server.

        :param file_type: The type of file to be downloaded.
        :type file_type: _TYPES_FILE ["gff", "gtf", "fasta"]
        :param sequence_nature: The nature of the sequence.
        :type sequence_nature: _TYPES_FILE_SEQ["dna", "ncrna"]
        :return: A tuple containing the path to the downloaded file, the annotation release version, and the assembly name.
        :rtype: Tuple[str, str, str]
        """
        self._check_file_type(file_type)
        self._check_sequence_nature_type(sequence_nature)

        ftp_directory, ftp_file = self._get_params(file_type, sequence_nature)
        dowloaded_file = self._download_and_decompress(self.ftp_link, ftp_directory, ftp_file)

        self.assembly_name = re.search("\\.([^\\.]*)\\.", Path(dowloaded_file).name).group().replace(".", "")

        return dowloaded_file, self.annotation_release, self.assembly_name

    def _get_params(self, file_type: _TYPES_FILE, sequence_nature: _TYPES_FILE_SEQ) -> Tuple[str, str]:
        """
        Constructs the FTP directory path and file name based on the file type and sequence nature.

        :param file_type: The type of file to be downloaded.
        :type file_type: _TYPES_FILE ["gff", "gtf", "fasta"]
        :param sequence_nature: The nature of the sequence.
        :type sequence_nature: _TYPES_FILE_SEQ["dna", "ncrna"]
        :return: A tuple containing the FTP directory path and the file name.
        :rtype: Tuple[str, str]
        """
        Path(self.dir_output).mkdir(parents=True, exist_ok=True)

        if self.annotation_release == "current":
            file_readme = self._download(self.ftp_link, "pub/", "current_README")
            with open(file_readme, "r") as handle:
                for line in handle:
                    if line.startswith("The current release is"):
                        self.annotation_release = line.strip().split("Ensembl ")[1]
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
    """
    A class for downloading and processing genomic data files from the NCBI FTP server.

    The `FtpLoaderNCBI` class handles the retrieval of genomic data files, such as GFF, GTF, and FASTA, from the NCBI FTP server.
    It supports downloading, decompressing, and mapping chromosome names based on the specified taxon, species, and annotation release.
    The class also manages the correct retrieval paths and handles different versions and structures of NCBI directories.

    :param taxon: The taxonomic group of the species (e.g., "vertebrate_mammalian").
    :type taxon: str
    :param species: The species name (e.g., "homo_sapiens").
    :type species: str
    :param annotation_release: The annotation release version to download (e.g., "109" or "current").
    :type annotation_release: str
    """

    def __init__(self, dir_output: str, taxon: str, species: str, annotation_release: str) -> None:
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

    def download_files(self, file_type: _TYPES_FILE) -> Tuple[str, str, str]:
        """
        Downloads the specified file type from the NCBI FTP server, decompresses it, and applies necessary chromosome name mappings.

        :param file_type: The type of file to be downloaded.
        :type file_type: _TYPES_FILE ["gff", "gtf", "fasta"]
        :return: A tuple containing the path to the downloaded file, the annotation release version, and the assembly name.
        :rtype: Tuple[str, str, str]
        """
        self._check_file_type(file_type)

        ftp_directory, ftp_file, ftp_file_chr_mapping = self._get_params(file_type)

        mapping = self._download_mapping_chr_names(ftp_directory, ftp_file_chr_mapping)
        dowloaded_file = self._download_and_decompress(self.ftp_link, ftp_directory, ftp_file)

        self.file_type_function[file_type](dowloaded_file, mapping)

        return dowloaded_file, self.annotation_release, self.assembly_name

    def _get_params(self, file_type: _TYPES_FILE) -> Tuple[str, str, str]:
        """
        Generates the necessary FTP directory paths and file names for downloading files from NCBI.

        :param file_type: The type of file to be downloaded.
        :type file_type: _TYPES_FILE ["gff", "gtf", "fasta"]
        :return: A tuple containing the FTP directory path, the file name, and the chromosome name mapping file name.
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

    def _download_mapping_chr_names(self, ftp_directory: str, ftp_file_chr_mapping: str) -> dict:
        """
        Downloads and processes the chromosome name mapping file from the NCBI FTP server.

        This function retrieves the chromosome name mapping file, skips comment lines, and extracts the mapping of RefSeq accession numbers to chromosome or scaffold names.
        It returns a dictionary where keys are RefSeq accession numbers, and values are the corresponding chromosome or scaffold names.

        :param ftp_directory: The FTP directory path where the mapping file is located.
        :type ftp_directory: str
        :param ftp_file_chr_mapping: The filename of the chromosome mapping file on the FTP server.
        :type ftp_file_chr_mapping: str
        :return: A dictionary mapping RefSeq accession numbers to chromosome or scaffold names.
        :rtype: dict
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

    def _map_chr_names_gene_annotation(self, ftp_file: str, mapping: dict) -> None:
        """
        Maps chromosome names in a gene annotation file from RefSeq accession numbers to chromosome names.

        This function reads a gene annotation file, replaces the sequence identifiers (seqid) with mapped chromosome names from a provided dictionary,
        and writes the updated annotations back to the file.

        :param ftp_file: The path to the gene annotation file (GTF/GFF).
        :type ftp_file: str
        :param mapping: A dictionary mapping RefSeq accession numbers to chromosome names.
        :type mapping: dict
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

    def _map_chr_names_genome_sequence(self, ftp_file: str, mapping: dict) -> None:
        """
        Maps chromosome names in a genome sequence file from RefSeq accession numbers to chromosome names.

        This function reads a genome sequence file, replaces sequence identifiers with mapped chromosome names from a provided dictionary,
        and writes the updated sequences back to the file.

        :param ftp_file: The path to the genome sequence file (FASTA).
        :type ftp_file: str
        :param mapping: A dictionary mapping RefSeq accession numbers to chromosome names.
        :type mapping: dict
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
