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
from typing import IO, get_args

import pandas as pd
from Bio import SeqIO

from oligo_designer_toolsuite._constants import _TYPES_FILE, _TYPES_FILE_SEQ
from oligo_designer_toolsuite._exceptions import ConfigurationError
from oligo_designer_toolsuite.utils import logger

############################################
# FTP Classes
############################################


class BaseFtpLoader:
    """
    A base class for downloading files via FTP and postprocessing the downloaded files.

    :param dir_output: Directory path where output files will be saved.
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

        if file_output is None:
            raise FileNotFoundError(f"File '{file_name}' not found on FTP server.")

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

    :param dir_output: Directory path where output files will be saved.
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
        self.assembly_name = ""
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
    ) -> tuple[str, str, str]:
        """
        Downloads and decompresses genomic data files from the Ensembl FTP server.

        :param file_type: The type of file to be downloaded.
        :type file_type: _TYPES_FILE ["gff", "gtf", "fasta"]
        :param sequence_nature: The nature of the sequence.
        :type sequence_nature: _TYPES_FILE_SEQ["dna", "ncrna"]
        :return: A tuple containing the path to the downloaded file, the annotation release version, and the assembly name.
        :rtype: tuple[str, str, str]
        """
        self._check_file_type(file_type)
        self._check_sequence_nature_type(sequence_nature)

        ftp_directory, ftp_file = self._get_params(file_type, sequence_nature)
        dowloaded_file = self._download_and_decompress(self.ftp_link, ftp_directory, ftp_file)

        match = re.search("\\.([^\\.]*)\\.", Path(dowloaded_file).name)
        if match is not None:
            self.assembly_name = match.group().replace(".", "")

        return dowloaded_file, self.annotation_release, self.assembly_name

    def _get_params(self, file_type: _TYPES_FILE, sequence_nature: _TYPES_FILE_SEQ) -> tuple[str, str]:
        """
        Constructs the FTP directory path and file name based on the file type and sequence nature.

        :param file_type: The type of file to be downloaded.
        :type file_type: _TYPES_FILE ["gff", "gtf", "fasta"]
        :param sequence_nature: The nature of the sequence.
        :type sequence_nature: _TYPES_FILE_SEQ["dna", "ncrna"]
        :return: A tuple containing the FTP directory path and the file name.
        :rtype: tuple[str, str]
        """
        Path(self.dir_output).mkdir(parents=True, exist_ok=True)

        if self.annotation_release == "current":
            file_version = self._download(self.ftp_link, "pub/", "VERSION")
            with open(file_version, "r") as handle:
                self.annotation_release = handle.readline().strip()
            os.remove(file_version)

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

    :param mode: NCBI parameter mode. Supported values are "species" and "assembly".
    :type mode: str | None
    :param taxon: The taxonomic group of the species (e.g., "vertebrate_mammalian").
    :type taxon: str | None
    :param species: The species name (e.g., "homo_sapiens").
    :type species: str | None
    :param annotation_release: The annotation release version to download (e.g., "109" or "current").
    :type annotation_release: str | None
    :param assembly_source: NCBI source selection. Supported values are "auto", "annotation_releases",
        "latest_assembly_versions", and "reference".
    :type assembly_source: str
    :param refseq_assembly_accession: Optional direct RefSeq assembly accession (e.g., "GCF_000001405.38").
    :type refseq_assembly_accession: str | None
    :param assembly_name: Optional direct assembly name (e.g., "GRCh38.p12").
    :type assembly_name: str | None
    """

    SUPPORTED_TAXA_SOURCES: dict[str, set[str]] = {
        "archaea": {"latest_assembly_versions", "reference"},
        "bacteria": {"latest_assembly_versions", "reference"},
        "fungi": {"latest_assembly_versions", "reference"},
        "invertebrate": {"annotation_releases", "latest_assembly_versions", "reference"},
        "metagenomes": {"latest_assembly_versions"},
        "plant": {"annotation_releases", "latest_assembly_versions", "reference"},
        "plants": {"annotation_releases", "latest_assembly_versions", "reference"},
        "protozoa": {"latest_assembly_versions", "reference"},
        "unknown": {"latest_assembly_versions"},
        "vertebrate_mammalian": {"annotation_releases", "latest_assembly_versions", "reference"},
        "vertebrate_other": {"annotation_releases", "latest_assembly_versions", "reference"},
        "viral": {"latest_assembly_versions"},
    }
    UNSUPPORTED_TAXA: set[str] = {"mitochondrion", "plasmids", "plastid"}
    # Determines how the assembly for a species is selected from the possible sources within the NCBI FTP directory.
    # 'annotation_releases' directory, should exist for all eukaryotic species and contains assemblies annotated with different annotation versions and the annotation version can be specified by 'annotation_release'.
    # 'latest_assembly_versions' directory is available for all species and contains the latest assembly.
    # 'reference' directory contains the reference genome. This is only available for a subset of species.
    # 'auto' automatically selects an assembly source in the following order (if available): 'annotation_releases', 'latest_assembly_versions'
    ALLOWED_ASSEMBLY_SOURCES: set[str] = {
        "auto",
        "annotation_releases",
        "latest_assembly_versions",
        "reference",
    }
    ALLOWED_MODES: set[str] = {"species", "assembly"}

    def __init__(
        self,
        dir_output: str,
        mode: str | None = None,
        taxon: str | None = None,
        species: str | None = None,
        annotation_release: str | None = None,
        assembly_source: str = "auto",
        refseq_assembly_accession: str | None = None,
        assembly_name: str | None = None,
    ) -> None:
        """Constructor for the FtpLoaderNCBI class."""
        super().__init__(dir_output)
        self.mode = mode
        self.taxon = taxon
        self.species = species
        self.annotation_release = annotation_release
        self.annotation_release_from_ncbi = "no_annotation_info"
        self.assembly_source = assembly_source
        self.assembly_name = assembly_name if assembly_name is not None else ""
        self.assembly_accession = refseq_assembly_accession if refseq_assembly_accession is not None else ""

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
        self._validate_mode_and_normalize_params()
        self._validate_taxon_and_source()

    def download_files(self, file_type: _TYPES_FILE) -> tuple[str, str, str]:
        """
        Downloads the specified file type from the NCBI FTP server, decompresses it, and applies necessary chromosome name mappings.

        :param file_type: The type of file to be downloaded.
        :type file_type: _TYPES_FILE ["gff", "gtf", "fasta"]
        :return: A tuple containing the path to the downloaded file, the annotation release version (as determined from information provided by NCBI), and the assembly name.
        :rtype: tuple[str, str, str]
        """
        self._check_file_type(file_type)

        ftp_directory, ftp_file, ftp_file_chr_mapping = self._get_params(file_type)

        mapping = self._download_mapping_chr_names(ftp_directory, ftp_file_chr_mapping)
        dowloaded_file = self._download_and_decompress(self.ftp_link, ftp_directory, ftp_file)

        self.file_type_function[file_type](dowloaded_file, mapping)

        return dowloaded_file, self.annotation_release_from_ncbi, self.assembly_name

    def _get_params(self, file_type: _TYPES_FILE) -> tuple[str, str, str]:
        """
        Generates the necessary FTP directory paths and file names for downloading files from NCBI.

        :param file_type: The type of file to be downloaded.
        :type file_type: _TYPES_FILE ["gff", "gtf", "fasta"]
        :return: A tuple containing the FTP directory path, the file name, and the chromosome name mapping file name.
        :rtype: tuple[str, str, str]
        """
        Path(self.dir_output).mkdir(parents=True, exist_ok=True)

        if self.mode == "assembly":
            ftp_directory = self._resolve_directory_from_direct_assembly()
            self._resolve_annotation_metadata(ftp_directory)
        else:
            source_subdir = self._resolve_source_subdir()
            ftp_directory = self._resolve_base_directory(source_subdir)
            self._resolve_assembly_metadata(ftp_directory)
            # all the numeric annotations (up until 110) have the following directory structure:
            # /genomes/refseq/vertebrate_mammalian/Homo_sapiens/annotation_releases/110
            # |- GCF_000001405.40_GRCh38.p14/
            # |- GCF_009914755.1_T2T-CHM13v2.0/
            # |- Homo_sapiens_AR110_annotation_report.xml
            # |- README_Homo_sapiens_annotation_release_110
            # there the assembly accession & name was inferred from README_Homo_sapiens_annotation_release_110
            # and one need to change into the GCF_000001405.40_GRCh38.p14/ directory where all the files are
            # annotations that are not numeric have the following directory structure:
            # /genomes/refseq/vertebrate_mammalian/Homo_sapiens/annotation_releases/GCF_009914755.1-RS_2025_08
            # |- GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz
            # |- ...
            # ]- README_GCF_009914755.1-RS_2025_08
            # there the assembly accession & name was inferred from README_GCF_009914755.1-RS_2025_08 directly in
            # this directory and we don't need to change it anymore.
            ftp_directory = self._resolve_assembly_directory(ftp_directory)
            self._resolve_annotation_metadata(ftp_directory)

        ftp_file = f"{self.assembly_accession}_{self.assembly_name}_{self.file_type_ending[file_type]}"
        ftp_file_chr_mapping = f"{self.assembly_accession}_{self.assembly_name}_assembly_report.txt"

        return ftp_directory, ftp_file, ftp_file_chr_mapping

    def _validate_mode_and_normalize_params(self) -> None:
        if self.mode not in self.ALLOWED_MODES:
            allowed_modes = ", ".join(sorted(self.ALLOWED_MODES))
            raise ConfigurationError(
                f"mode '{self.mode}' is not supported. Supported values are: {allowed_modes}."
            )

        if self.mode == "assembly":
            has_direct_accession = bool(self.assembly_accession)
            has_direct_name = bool(self.assembly_name)
            if not has_direct_accession or not has_direct_name:
                raise ConfigurationError(
                    "In mode='assembly', both 'refseq_assembly_accession' and 'assembly_name' must be provided."
                )
            if any(param is not None for param in (self.taxon, self.species, self.annotation_release)):
                raise ConfigurationError(
                    "In mode='assembly', taxon/species/annotation_release must not be provided."
                )
            if self.assembly_source != "auto":
                raise ConfigurationError(
                    "In mode='assembly', assembly_source cannot be set and must remain 'auto'."
                )
            self.annotation_release = "unknown"
        elif self.mode == "species":
            if bool(self.assembly_accession) or bool(self.assembly_name):
                raise ConfigurationError(
                    "In mode='species', refseq_assembly_accession/assembly_name must not be provided."
                )
            if self.taxon is None or self.species is None or self.annotation_release is None:
                raise ConfigurationError(
                    "In mode='species', 'taxon', 'species', and 'annotation_release' must be provided."
                )

    def _validate_taxon_and_source(self) -> None:
        if self.mode == "assembly":
            return

        if self.taxon is None:
            raise ConfigurationError("Taxon is not defined.")
        if self.taxon in self.UNSUPPORTED_TAXA:
            raise ConfigurationError(f"Taxon '{self.taxon}' is not supported for automated NCBI downloads.")
        if self.taxon not in self.SUPPORTED_TAXA_SOURCES:
            supported_taxa = ", ".join(sorted(self.SUPPORTED_TAXA_SOURCES))
            raise ConfigurationError(
                f"Taxon '{self.taxon}' is not supported. Supported taxa are: {supported_taxa}."
            )
        if self.assembly_source not in self.ALLOWED_ASSEMBLY_SOURCES:
            allowed_sources = ", ".join(sorted(self.ALLOWED_ASSEMBLY_SOURCES))
            raise ConfigurationError(
                f"assembly_source '{self.assembly_source}' is not supported. "
                f"Supported values are: {allowed_sources}."
            )

        available_sources = self.SUPPORTED_TAXA_SOURCES[self.taxon]
        if self.assembly_source != "auto" and self.assembly_source not in available_sources:
            supported_sources = ", ".join(sorted(available_sources))
            raise ConfigurationError(
                f"assembly_source '{self.assembly_source}' is not available for taxon '{self.taxon}'. "
                f"Supported sources for this taxon are: {supported_sources}."
            )

        if self.annotation_release != "current" and (
            self.assembly_source in {"latest_assembly_versions", "reference"}
            or (self.assembly_source == "auto" and "annotation_releases" not in available_sources)
        ):
            raise ConfigurationError(
                "A numeric annotation_release is only supported with assembly_source='annotation_releases'. "
                "Use annotation_release='current' for latest_assembly_versions/reference."
            )

    def _resolve_source_subdir(self) -> str:
        if self.taxon is None:
            raise ConfigurationError("taxon cannot be none.")
        available_sources = self.SUPPORTED_TAXA_SOURCES[self.taxon]
        if self.assembly_source != "auto":
            return self.assembly_source
        if "annotation_releases" in available_sources:
            return "annotation_releases"
        if "latest_assembly_versions" in available_sources:
            return "latest_assembly_versions"
        supported_sources = ", ".join(sorted(available_sources))
        raise ConfigurationError(
            f"Taxon '{self.taxon}' has no supported automatic source. Available sources: {supported_sources}."
        )

    def _resolve_base_directory(self, source_subdir: str) -> str:
        if self.taxon is None or self.species is None:
            raise ConfigurationError("Taxon/species are not defined.")
        base_directory = f"genomes/refseq/{self.taxon}/{self.species}/{source_subdir}/"
        if source_subdir != "annotation_releases":
            if self.annotation_release != "current":
                raise ConfigurationError(
                    "annotation_release must be 'current' when using assembly_source " f"'{source_subdir}'."
                )
            entries = self._list_ftp_entries(base_directory)
            gcf_entry = next((entry for entry in entries if entry.startswith("GCF")), None)
            if gcf_entry is None:
                if source_subdir != "reference":
                    additional_error_msg = " Try 'reference' for 'assembly_source'."
                else:
                    additional_error_msg = ""
                raise ConfigurationError(
                    f"No GCF assembly folder found in '{base_directory}' on NCBI FTP.{additional_error_msg}"
                )
            return f"{base_directory}{gcf_entry}/"

        if self.annotation_release == "current":
            current_directory = f"{base_directory}current/"
            # inside current dir there is the directory containing the data for the lateest annotation release
            # in case there are several directories, use the first one
            entries = self._list_ftp_entries(current_directory)
            if not entries:
                raise ConfigurationError(
                    f"Could not resolve current annotation release in '{current_directory}'."
                )
            self.annotation_release = entries[0]
        return f"{base_directory}{self.annotation_release}/"

    def _resolve_directory_from_direct_assembly(self) -> str:
        if not self.assembly_accession or not self.assembly_name:
            raise ConfigurationError(
                "Direct assembly mode requires both 'refseq_assembly_accession' and 'assembly_name'."
            )
        match = re.match(r"^(GCF)_(\d+)\.(\d+)$", self.assembly_accession)
        if match is None:
            raise ConfigurationError(
                f"Invalid RefSeq assembly accession '{self.assembly_accession}'. "
                "Expected format like 'GCF_000001405.38'."
            )
        accession_prefix = match.group(1)
        accession_number = match.group(2)
        if len(accession_number) != 9:
            raise ConfigurationError(
                f"Invalid RefSeq assembly accession '{self.assembly_accession}'. "
                "Expected format like 'GCF_000001405.38'."
            )
        chunks = [accession_number[i : i + 3] for i in range(0, len(accession_number), 3)]
        grouped_number_path = "/".join(chunks)
        return (
            f"genomes/all/{accession_prefix}/{grouped_number_path}/"
            f"{self.assembly_accession}_{self.assembly_name}/"
        )

    def _resolve_assembly_metadata(self, ftp_directory: str) -> None:
        self.assembly_name = ""
        self.assembly_accession = ""

        # all species where information about the used annotation is available
        # have a README_{annotation name} file; either in the directory that contains
        # the subdirectory with the data and also in the directory with all the data
        # define parser helpers
        def parse_readme(handle: IO) -> None:
            for line in handle:
                if line.startswith("ASSEMBLY NAME:"):
                    self.assembly_name = line.strip().split("\t", 1)[1]
                if line.startswith("ASSEMBLY ACCESSION:"):
                    self.assembly_accession = line.strip().split("\t", 1)[1]
                    break

        def parse_assembly_report(handle: IO) -> None:
            for line in handle:
                if line.startswith("# Assembly name:"):
                    self.assembly_name = line.split(":", 1)[1].strip()
                if line.startswith("# RefSeq assembly accession:"):
                    self.assembly_accession = line.split(":", 1)[1].strip()
                    break

        candidates = [
            (r"README_(?!patch_release\.txt$).*", parse_readme),
            (r".*_assembly_report\.txt", parse_assembly_report),
        ]
        success = False
        for pattern, parser in candidates:
            try:
                file_path = self._download(self.ftp_link, ftp_directory, pattern)
                try:
                    with open(file_path, "r") as handle:
                        parser(handle)
                finally:
                    os.remove(file_path)
                if self.assembly_accession and self.assembly_name:
                    success = True
                    break
            except FileNotFoundError:
                # try next candidate
                continue

        if not success:
            raise ConfigurationError(
                f"Could not resolve assembly metadata from '{ftp_directory}' on NCBI FTP."
            )

    def _resolve_annotation_metadata(self, ftp_directory: str) -> None:
        """
        All species annotated with the NCBI Eukaryotic Genome Annotation Pipeline should have a
        README_{annotation name} file that contains the annotation name used by NCBI. In case
        this is not found, keep the default value ("no_annotation_info").
        """
        try:
            file_readme = self._download(self.ftp_link, ftp_directory, r"README_(?!patch_release\.txt$).*")
            with open(file_readme, "r") as handle:
                for line in handle:
                    if line.startswith("ANNOTATION RELEASE NAME:"):
                        annotation_release_name = line.split(":", 1)[1].strip()
                        self.annotation_release_from_ncbi = annotation_release_name.replace(" ", "_")
                        break
            os.remove(file_readme)
        except FileNotFoundError:
            logger.warning("No annotation name information available from NCBI.")

    def _resolve_assembly_directory(self, ftp_directory: str) -> str:
        assembly_dir = f"{ftp_directory}{self.assembly_accession}_{self.assembly_name}"
        assembly_dir = assembly_dir.replace(" ", "_")
        ftp = FTP(self.ftp_link)
        ftp.login()
        try:
            ftp.cwd(assembly_dir)
            return assembly_dir
        except error_perm:
            return ftp_directory
        finally:
            ftp.quit()

    def _list_ftp_entries(self, ftp_directory: str) -> list[str]:
        ftp = FTP(self.ftp_link)
        ftp.login()
        ftp.cwd(ftp_directory)
        entries = ftp.nlst()
        ftp.quit()
        return entries

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
            names_list = names[1:].split()

        assembly_report = pd.read_table(file_mapping, names=names_list, sep="\t", comment="#")

        mapping_chromosome_df = assembly_report[assembly_report["Sequence-Role"] == "assembled-molecule"]
        mapping_chromosome = pd.Series(
            mapping_chromosome_df["Sequence-Name"].values,
            index=mapping_chromosome_df["RefSeq-Accn"],
        ).to_dict()

        mapping_scaffolds_df = assembly_report[assembly_report["Sequence-Role"] != "assembled-molecule"]
        mapping_scaffolds = pd.Series(
            mapping_scaffolds_df["GenBank-Accn"].values,
            index=mapping_scaffolds_df["RefSeq-Accn"],
        ).to_dict()

        mapping: dict[str, str] = mapping_chromosome
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
                    logger.warning("No mapping for accession number: {}".format(accession_number))

        os.replace(file_tmp, ftp_file)
