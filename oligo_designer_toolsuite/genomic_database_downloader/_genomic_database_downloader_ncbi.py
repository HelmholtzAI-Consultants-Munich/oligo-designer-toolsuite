import hashlib
import itertools
import os
import re
from abc import ABC, abstractmethod
from dataclasses import dataclass
from ftplib import FTP, error_perm
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from typing_extensions import IO, ClassVar, override

from oligo_designer_toolsuite._constants import _TYPES_FILE, _TYPES_FILE_SEQ
from oligo_designer_toolsuite._exceptions import ConfigurationError
from oligo_designer_toolsuite.genomic_database_downloader._genomic_database_downloader_base import (
    FetchedGenomicEntity,
    FtpUtil,
    GenomicDatabaseConfigBase,
    GenomicDatabaseDownloaderBase,
    GenomicEntityContext,
    RemoteInfo,
)
from oligo_designer_toolsuite.utils import logger


@dataclass(frozen=True)
class NcbiGenomicEntity:
    mode: str | None
    taxon: str | None
    species: str | None
    annotation_release: str | None
    assembly_name: str | None
    assembly_accession: str | None
    assembly_source: str | None


class NcbiDownloadMode(ABC):
    """
    Abstract interface for resolving NCBI assemblies and FTP locations.
    """

    @abstractmethod
    def get_ftp_directory(self) -> tuple[str, str, str]:
        """
        Resolve the FTP directory and assembly metadata.

        :return: FTP directory, assembly accession and assembly name.
        :rtype: tuple[str, str, str]
        """

    @abstractmethod
    def validate_and_set(self, entity: NcbiGenomicEntity) -> None:
        """
        Validate an entity configuration and initialize strategy state.

        :param entity: Entity configuration to validate.
        :type entity: NcbiGenomicEntity
        """

    @abstractmethod
    def _validate_taxon_and_source(self) -> None:
        """
        Validate taxon and assembly source configuration.

        :raises ConfigurationError: If the configuration is invalid.
        """


class AssemblyMode(NcbiDownloadMode):
    """
    Resolution strategy for direct assembly-based downloads.

    This strategy resolves FTP locations directly from a RefSeq assembly
    accession and assembly name.

    :param entity: Entity describing the target assembly.
    :type entity: NcbiGenomicEntity
    """

    def __init__(self, entity: NcbiGenomicEntity):
        """Constructor for AssemblyMode class."""
        self.validate_and_set(entity)

    def _validate_taxon_and_source(self) -> None:
        """
        Not Implemented, since Assembly Mode does not need the taxon or assembly source.
        """

    def _resolve_directory_from_direct_assembly(self) -> str:
        """
        Build the FTP directory for a RefSeq assembly.

        :return: FTP directory containing assembly resources.
        :rtype: str

        :raises ConfigurationError: If the assembly accession is invalid.
        """
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

    def get_ftp_directory(self) -> tuple[str, str, str]:
        """
        Resolve FTP directory information for the configured assembly and adds assembly accession and assembly name for easier retrieval.

        :return: FTP directory, assembly accession and assembly name.
        :rtype: tuple[str, str, str]
        """
        return (
            self._resolve_directory_from_direct_assembly(),
            self.assembly_accession,
            self.assembly_name,
        )

    def validate_and_set(self, entity: NcbiGenomicEntity) -> None:
        """
        Validate and initialize assembly mode parameters.

        :param entity: Entity configuration.
        :type entity: NcbiGenomicEntity

        :raises ConfigurationError: If assembly mode parameters are invalid.
        """
        if entity.assembly_accession is None or entity.assembly_name is None:
            raise ConfigurationError(
                "In mode='assembly', both 'refseq_assembly_accession' and 'assembly_name' must be provided."
            )
        if any(param is not None for param in (entity.taxon, entity.species, entity.annotation_release)):
            raise ConfigurationError(
                "In mode='assembly', taxon/species/annotation_release must not be provided."
            )
        if entity.assembly_source != "auto":
            if entity.assembly_source is None:
                self.assembly_source = "auto"
            else:
                raise ConfigurationError(
                    "In mode='assembly', assembly_source cannot be set and must remain 'auto'."
                )

        self.annotation_release = "unknown"
        self.assembly_accession = entity.assembly_accession
        self.assembly_name = entity.assembly_name


class SpeciesMode(NcbiDownloadMode):
    """
    Resolution strategy for species-based downloads.

    This strategy resolves species, annotation releases and assembly
    information through the NCBI RefSeq FTP structure.

    :param host: FTP host name.
    :type host: str
    :param entity: Entity configuration.
    :type entity: NcbiGenomicEntity
    :param ftp_util: FTP helper instance.
    :type ftp_util: FtpUtil
    """

    SUPPORTED_TAXA_SOURCES: ClassVar[dict[str, set[str]]] = {
        "archaea": {"latest_assembly_versions", "reference"},
        "bacteria": {"latest_assembly_versions", "reference"},
        "fungi": {"latest_assembly_versions", "reference"},
        "invertebrate": {
            "annotation_releases",
            "latest_assembly_versions",
            "reference",
        },
        "metagenomes": {"latest_assembly_versions"},
        "plant": {"annotation_releases", "latest_assembly_versions", "reference"},
        "plants": {"annotation_releases", "latest_assembly_versions", "reference"},
        "protozoa": {"latest_assembly_versions", "reference"},
        "unknown": {"latest_assembly_versions"},
        "vertebrate_mammalian": {
            "annotation_releases",
            "latest_assembly_versions",
            "reference",
        },
        "vertebrate_other": {
            "annotation_releases",
            "latest_assembly_versions",
            "reference",
        },
        "viral": {"latest_assembly_versions"},
    }
    UNSUPPORTED_TAXA: ClassVar[set[str]] = {"mitochondrion", "plasmids", "plastid"}
    # Determines how the assembly for a species is selected from the possible sources within the NCBI FTP directory.
    # 'annotation_releases' directory, should exist for all eukaryotic species and contains assemblies annotated with different annotation versions and the annotation version can be specified by 'annotation_release'.
    # 'latest_assembly_version' directory is available for all species and contains the latest assembly.
    # 'reference' directory contains the reference genome. This is only available for a subset of species.
    # 'auto' automatically selects an assembly source in the following order (if available): 'annotation_releases', 'latest_assembly_version'
    ALLOWED_ASSEMBLY_SOURCES: ClassVar[set[str]] = {
        "auto",
        "annotation_releases",
        "latest_assembly_versions",
        "reference",
    }

    def __init__(self, host: str, entity: NcbiGenomicEntity, ftp_util: FtpUtil):
        """Constructor for AssemblyMode class."""
        self.validate_and_set(entity)

        self.host = host
        self.ftp_util = ftp_util

    def validate_and_set(self, entity: NcbiGenomicEntity) -> None:
        """
        Validate and initialize species mode parameters.

        :param entity: Entity configuration.
        :type entity: NcbiGenomicEntity

        :raises ConfigurationError: If required parameters are missing or invalid.
        """
        if bool(entity.assembly_accession) or bool(entity.assembly_name):
            raise ConfigurationError(
                "In mode='species', refseq_assembly_accession/assembly_name must not be provided."
            )

        if entity.taxon is None:
            raise ConfigurationError("No taxon defined.")

        if entity.species is None:
            raise ConfigurationError("No species defined.")

        if entity.annotation_release is None:
            annotation_release = "current"
            logger.warning(f"No annotation release defined. Using default release {annotation_release}!")

        if entity.assembly_source is None:
            self.assembly_source = "auto"
            logger.warning("No assembly source defined. Using default assembly source auto!")
        else:
            self.assembly_source = entity.assembly_source

        self.taxon = entity.taxon
        self.species = entity.species
        self.annotation_release = entity.annotation_release

    def get_ftp_directory(self) -> tuple[str, str, str]:
        """
        Resolve the FTP directory and assembly metadata.

        :return: FTP directory, assembly accession and assembly name.
        :rtype: tuple[str, str, str]
        """

        source_subdir = self._resolve_source_subdir()
        ftp_dir = self._resolve_base_directory(source_subdir)
        assembly_accession, assembly_name = self._resolve_assembly_metadata(ftp_dir)
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
        return (
            self._resolve_assembly_directory(ftp_dir, assembly_accession, assembly_name),
            assembly_accession,
            assembly_name,
        )

    def _resolve_assembly_directory(
        self, ftp_directory: str, assembly_accession: str, assembly_name: str
    ) -> str:
        """
        Resolve the directory containing assembly files.

        :param ftp_directory: Base FTP directory.
        :type ftp_directory: str
        :param assembly_accession: RefSeq assembly accession.
        :type assembly_accession: str
        :param assembly_name: Assembly name.
        :type assembly_name: str
        :return: Final FTP directory containing downloadable resources.
        :rtype: str
        """
        assembly_dir = f"{ftp_directory}{assembly_accession}_{assembly_name}"
        assembly_dir = assembly_dir.replace(" ", "_")
        ftp = FTP(self.host)
        ftp.login()
        try:
            ftp.cwd(assembly_dir)
            return assembly_dir
        except error_perm:
            return ftp_directory
        finally:
            ftp.quit()

    def _resolve_source_subdir(self) -> str:
        """
        Resolve the assembly source subdirectory.

        :return: Selected source subdirectory.
        :rtype: str
        :raises ConfigurationError: If no valid source can be determined.
        """
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
        """
        Resolve the FTP base directory for the configured entity.

        :param source_subdir: Selected assembly source.
        :type source_subdir: str
        :return: FTP directory containing assembly information.
        :rtype: str
        :raises ConfigurationError: If the directory configuration is invalid.
        """
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

    def _resolve_assembly_metadata(self, ftp_directory: str) -> tuple[str, str]:
        """
        Resolve assembly accession and assembly name from FTP metadata files.

        :param ftp_directory: FTP directory to inspect.
        :type ftp_directory: str
        :return: Assembly accession and assembly name.
        :rtype: tuple[str, str]
        :raises ConfigurationError: If assembly metadata cannot be determined.
        """
        assembly_accession = None
        assembly_name = None

        # all species where information about the used annotation is available
        # have a README_{annotation name} file; either in the directory that contains
        # the subdirectory with the data and also in the directory with all the data
        # define parser helpers
        def parse_readme(handle: IO) -> tuple[str | None, str | None]:
            assembly_accession = None
            assembly_name = None

            for line in handle:
                if line.startswith("ASSEMBLY NAME:"):
                    assembly_name = line.strip().split("\t", 1)[1]
                if line.startswith("ASSEMBLY ACCESSION:"):
                    assembly_accession = line.strip().split("\t", 1)[1]
                    break

            return assembly_accession, assembly_name

        def parse_assembly_report(handle: IO) -> tuple[str | None, str | None]:
            assembly_accession = None
            assembly_name = None

            for line in handle:
                if line.startswith("# Assembly name:"):
                    assembly_name = line.split(":", 1)[1].strip()
                if line.startswith("# RefSeq assembly accession:"):
                    assembly_accession = line.split(":", 1)[1].strip()
                    break

            return assembly_accession, assembly_name

        candidates = [
            (r"README_(?!patch_release\.txt$).*", parse_readme),
            (r".*_assembly_report\.txt", parse_assembly_report),
        ]
        for pattern, parser in candidates:
            try:
                file_path = self.ftp_util.download(
                    ftp_directory,
                    self.ftp_util.get_exact_filename(ftp_directory, pattern),
                )
                try:
                    with open(file_path, "r") as handle:
                        assembly_accession, assembly_name = parser(handle)
                finally:
                    os.remove(file_path)
                if assembly_accession and assembly_name:
                    break
            except FileNotFoundError:
                # try next candidate
                continue

        if assembly_accession is None or assembly_name is None:
            raise ConfigurationError(
                f"Could not resolve assembly metadata from '{ftp_directory}' on NCBI FTP."
            )

        return assembly_accession, assembly_name

    def _validate_taxon_and_source(self) -> None:
        """
        Validate given taxon and assembly source.

        :raises ConfigurationError: If taxon or assembly source are either invalid or unsupported.
        """
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

    def _list_ftp_entries(self, ftp_directory: str) -> list[str]:
        """
        List entries contained in an FTP directory.

        :param ftp_directory: FTP directory to inspect.
        :type ftp_directory: str
        :return: Directory entries.
        :rtype: list[str]
        """
        ftp = FTP(self.host)
        ftp.login()
        ftp.cwd(ftp_directory)
        entries = ftp.nlst()
        ftp.quit()
        return entries


class GenomicDatabaseConfigNcbi(GenomicDatabaseConfigBase):
    name: ClassVar[str] = "ncbi"
    host: ClassVar[str] = "ftp.ncbi.nlm.nih.gov"
    base_path: ClassVar[str] = "/genomes/refseq"
    checksums_filename: ClassVar[str] = "md5checksums.txt"


class GenomicDatabaseDownloaderNcbi(GenomicDatabaseDownloaderBase[NcbiGenomicEntity]):
    """
    A class for downloading and processing genomic data files from the NCBI FTP server.

    The `GenomicDatabaseDownloaderNcbi` class handles the retrieval of genomic data files, such as GFF, GTF, and FASTA, from the NCBI FTP server.
    It supports downloading, decompressing, and mapping chromosome names based on the specified taxon, species, and annotation release.
    The class also manages the correct retrieval paths and handles different versions and structures of NCBI directories.

    :param dir_output: Local cache directory for downloaded files.
    :type dir_output: str
    """

    name: ClassVar[str] = GenomicDatabaseConfigNcbi.name
    host: ClassVar[str] = GenomicDatabaseConfigNcbi.host
    base_path: ClassVar[str] = GenomicDatabaseConfigNcbi.base_path
    checksums_filename: ClassVar[str] = GenomicDatabaseConfigNcbi.checksums_filename

    ALLOWED_MODES: set[str] = {"species", "assembly"}

    def __init__(
        self,
        dir_output: str,
    ) -> None:
        """Constructor for GenomicDatabaseDownloaderNcbi class."""
        super().__init__(Path(dir_output))

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

        self.annotation_release_from_ncbi: str | None = None

    @override
    def _matches_checksum(self, file_path: Path, expected_checksum: str) -> bool:
        """
        Check whether a file matches the expected MD5 checksum.

        :param file_path: Path to the local file.
        :type file_path: pathlib.Path
        :param expected_checksum: Expected checksum value.
        :type expected_checksum: str
        :return: True if the checksum matches, otherwise False.
        :rtype: bool
        """
        with open(file_path, "rb") as f:
            digest = hashlib.file_digest(f, "md5")

        return digest.hexdigest() == expected_checksum

    def fetch_genomic_entity(
        self,
        entity: NcbiGenomicEntity,
        files: list[_TYPES_FILE] | None = None,
        sequence_nature: _TYPES_FILE_SEQ = "dna",
    ) -> FetchedGenomicEntity:
        """
        Fetch a genomic entity from cache or download it.

        :param entity: Genomic entity to fetch.
        :type entity: NcbiGenomicEntity
        :param files: Requested file types.
        :type files: list[_TYPES_FILE] | None
        :param sequence_nature: Requested sequence nature.
        :type sequence_nature: _TYPES_FILE_SEQ
        :return: Downloaded files together with resolved metadata.
        :rtype: FetchedGenomicEntity
        :raises RuntimeError: If no cache directory is configured.
        :raises RuntimeError: If a required resource cannot be located.
        :raises RuntimeError: If a download fails.
        :raises RuntimeError: If checksum verification fails.
        """
        self.mode = self._get_mode(entity)
        self.mode._validate_taxon_and_source()

        fetched_genomic_entity = super().fetch_genomic_entity(entity, files)

        for file_type, file_path in fetched_genomic_entity["fetched_files"].items():
            remote_dir, _, mapping_file = self._get_params(file_type)
            self._apply_mapping_chr_names(file_type, file_path, remote_dir, mapping_file)

        return fetched_genomic_entity

    @override
    def _get_entity_context(
        self, entity: NcbiGenomicEntity, files: list[_TYPES_FILE], sequence_nature: _TYPES_FILE_SEQ
    ) -> GenomicEntityContext:
        """
        Retrieve download context for a genomic entity.

        :param entity: Genomic entity to resolve.
        :type entity: NcbiGenomicEntity
        :param files: Requested file types.
        :type files: list[_TYPES_FILE]
        :param sequence_nature: Requested sequence nature.
        :type sequence_nature: _TYPES_FILE_SEQ
        :return: Download context for the requested entity.
        :rtype: GenomicEntityContext
        """
        remote_file_infos: dict[_TYPES_FILE, RemoteInfo] = {}
        for file_type in files:
            file_type_remote_dir, file_type_remote_filename, _ = self._get_params(file_type)
            remote_file_infos[file_type] = RemoteInfo(file_type_remote_dir, file_type_remote_filename)

        annotation_release = (
            self.annotation_release_from_ncbi
            if self.annotation_release_from_ncbi is not None
            else "no_annotation_info"
        )

        return GenomicEntityContext(
            remote_file_infos,
            annotation_release,
            self.assembly_name,
            self.assembly_accession,
        )

    def _get_mode(self, entity: NcbiGenomicEntity) -> NcbiDownloadMode:
        """
        Create the resolution mode for an entity.

        :param entity: Entity configuration.
        :type entity: NcbiGenomicEntity
        :return: Resolution mode instance.
        :rtype: NcbiDownloadMode
        :raises ConfigurationError: If the configured mode is unsupported.
        """
        match entity.mode:
            case "assembly":
                return AssemblyMode(entity)
            case "species":
                return SpeciesMode(self.host, entity, self.ftp_util)
            case _:
                allowed_modes = ", ".join(sorted(self.ALLOWED_MODES))

                if entity.mode is None:
                    raise ConfigurationError("For source='ncbi', parameter 'mode' must be provided.")

                raise ConfigurationError(
                    f"mode '{entity.mode}' is not supported. Supported values for source='ncbi' are: {allowed_modes}."
                )

    def _apply_mapping_chr_names(
        self,
        file_type: _TYPES_FILE,
        ftp_file: str,
        ftp_directory: str,
        ftp_mapping_file: str,
    ) -> None:
        """
        Apply chromosome name mapping to a downloaded file.

        :param file_type: File type being processed.
        :type file_type: _TYPES_FILE
        :param ftp_file: Local file path.
        :type ftp_file: str
        :param ftp_directory: FTP directory containing mapping resources.
        :type ftp_directory: str
        :param ftp_mapping_file: Mapping filename.
        :type ftp_mapping_file: str
        """
        mapping = self._download_mapping_chr_names(ftp_directory, ftp_mapping_file)

        self.file_type_function[file_type](ftp_file, mapping)

    def _get_params(self, file_type: _TYPES_FILE) -> tuple[str, str, str]:
        """
        Generates the necessary FTP directory paths and file names for downloading files from NCBI.

        :param file_type: The type of file to be downloaded.
        :type file_type: _TYPES_FILE ["gff", "gtf", "fasta"]
        :return: A tuple containing the FTP directory path, the file name, and the chromosome name mapping file name.
        :rtype: tuple[str, str, str]
        """
        ftp_directory, assembly_accession, assembly_name = self.mode.get_ftp_directory()

        self._resolve_annotation_metadata(ftp_directory)

        self.assembly_accession = assembly_accession
        self.assembly_name = assembly_name

        ftp_file = f"{assembly_accession}_{assembly_name}_{self.file_type_ending[file_type]}"
        ftp_file_chr_mapping = f"{assembly_accession}_{assembly_name}_assembly_report.txt"

        return ftp_directory, ftp_file, ftp_file_chr_mapping

    def _resolve_annotation_metadata(self, ftp_directory: str) -> None:
        """
        Resolve annotation release metadata from NCBI README files.

        The resolved annotation release name is stored in
        ``annotation_release_from_ncbi`` if available.

        :param ftp_directory: FTP directory containing annotation metadata.
        :type ftp_directory: str
        """
        try:
            file_readme = self.ftp_util.download(
                ftp_directory,
                self.ftp_util.get_exact_filename(ftp_directory, r"README_(?!patch_release\.txt$).*"),
            )
            with open(file_readme, "r") as handle:
                for line in handle:
                    if line.startswith("ANNOTATION RELEASE NAME:"):
                        annotation_release_name = line.split(":", 1)[1].strip()
                        self.annotation_release_from_ncbi = annotation_release_name.replace(" ", "_")
                        break
            os.remove(file_readme)
        except FileNotFoundError:
            logger.warning("No annotation name information available from NCBI.")

    def _download_mapping_chr_names(self, ftp_directory: str, ftp_file_chr_mapping: str) -> dict:
        """Literal["fasta", "gtf", "gff"]
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
        file_mapping = self.ftp_util.download(ftp_directory, ftp_file_chr_mapping)

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
        file_tmp = os.path.join(str(self.cache_dir), "temp.gtf")

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
        file_tmp = os.path.join(str(self.cache_dir), "temp.fna")

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
