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

from oligo_designer_toolsuite._constants import _TYPES_FILE
from oligo_designer_toolsuite._exceptions import ConfigurationError
from oligo_designer_toolsuite.ftp_loader._ftp_loader_base import (
    FtpLoaderBase,
    FtpUtil,
    GenomicEntityContext,
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


class ModeStrategy(ABC):
    @abstractmethod
    def get_ftp_directory(self) -> tuple[str, str, str]:
        pass

    @abstractmethod
    def validate_and_set(self, entity: NcbiGenomicEntity) -> None:
        pass

    @abstractmethod
    def _validate_taxon_and_source(self) -> None:
        pass


class AssemblyMode(ModeStrategy):
    def __init__(self, entity: NcbiGenomicEntity):
        self.validate_and_set(entity)

    def _validate_taxon_and_source(self) -> None:
        pass

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

    def get_ftp_directory(self) -> tuple[str, str, str]:
        return (
            self._resolve_directory_from_direct_assembly(),
            self.assembly_accession,
            self.assembly_name,
        )

    def validate_and_set(self, entity: NcbiGenomicEntity) -> None:
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


class SpeciesMode(ModeStrategy):
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
        self.validate_and_set(entity)

        self.host = host
        self.ftp_util = ftp_util

    def validate_and_set(self, entity: NcbiGenomicEntity) -> None:
        if bool(entity.assembly_accession) or bool(entity.assembly_name):
            raise ConfigurationError(
                "In mode='species', refseq_assembly_accession/assembly_name must not be provided."
            )

        if entity.taxon is None:
            raise ConfigurationError(f"No taxon defined.")

        if entity.species is None:
            raise ConfigurationError(f"No species defined.")

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
        source_subdir = self._resolve_source_subdir()
        ftp_dir = self._resolve_base_directory(source_subdir)
        self._resolve_assembly_metadata(ftp_dir)
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
            self._resolve_assembly_directory(ftp_dir),
            self.assembly_accession,
            self.assembly_name,
        )

    def _resolve_assembly_directory(self, ftp_directory: str) -> str:
        assembly_dir = f"{ftp_directory}{self.assembly_accession}_{self.assembly_name}"
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

    def _resolve_assembly_metadata(self, ftp_directory: str) -> None:
        pass

        # all species where information about the used annotation is available
        # have a README_{annotation name} file; either in the directory that contains
        # the subdirectory with the data and also in the directory with all the data
        # define parser helpers
        def parse_readme(handle: IO) -> None:
            for line in handle:
                if line.startswith("ASSEMBLY NAME:"):
                    assembly_name = line.strip().split("\t", 1)[1]
                if line.startswith("ASSEMBLY ACCESSION:"):
                    assembly_accession = line.strip().split("\t", 1)[1]
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
                file_path = self.ftp_util.download(
                    ftp_directory,
                    self.ftp_util.get_exact_filename(ftp_directory, pattern),
                )
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

    def _validate_taxon_and_source(self) -> None:
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
        ftp = FTP(self.host)
        ftp.login()
        ftp.cwd(ftp_directory)
        entries = ftp.nlst()
        ftp.quit()
        return entries


class FtpLoaderNCBI(FtpLoaderBase[NcbiGenomicEntity]):
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

    """A genomic database to interface with Ncbi."""

    name: ClassVar[str] = "ncbi"
    host: ClassVar[str] = "ftp.ncbi.nlm.nih.gov"
    base_path: ClassVar[str] = "/genomes/refseq"
    checksums_filename: ClassVar[str] = "md5checksums.txt"

    ALLOWED_MODES: set[str] = {"species", "assembly"}

    # ---- Directory Discovery ----
    @override
    def _build_species_mapping(self, species_dirs: dict[str, list[str]]) -> dict[str, list[str]]:
        """Builds the result of fetch_species_mapping.

        Arguments:
            species_dirs {dict[str, list[str]]} -- dict mapping taxon directories to their respective species.

        Returns:
            dict[str, list[str]] -- dict mapping taxons to their respective species.
        """
        return species_dirs

    def _get_all_releases_dir(self, ftp: FTP, taxon: str, species: str) -> str | None:
        """Searches for the remote directory containing all releases for a specific species.

        Arguments:
            ftp {ftplib.FTP} -- Active FTP handler.
            taxon {str} -- Taxon of the species.
            species {str} -- Species to look up, must be contained in taxon.

        Returns:
            str | None -- Directory containing all releases or None if not found.
        """
        candidate_subdirs = ["annotation_releases", "all_assembly_versions"]
        species_base_dir = f"{self.base_path}/{taxon}/{species}"

        # Look up all subdirectories for the species
        ftp.cwd(species_base_dir)
        subdirs = set(self._get_subdirs(ftp))

        # Return first releases directory found
        for candidate_subdir in candidate_subdirs:
            if candidate_subdir in subdirs:
                return f"{species_base_dir}/{candidate_subdir}"
        return None

    def fetch_annotations_releases(self, taxon: str, species: str) -> list[str] | None:
        """Fetches all available annotation releases for a specific species.

        Arguments:
            taxon {str} -- Taxon of the species.
            species {str} -- Species, must be contained in taxon.

        Notes:
            This function is decorated with our generic cache to reduce refetching.

        Returns:
            list[str] | None -- list of available annotation releases or None if none were found.
        """

        with FTP(self.host) as ftp:
            ftp.login()
            release_dir = self._get_all_releases_dir(ftp, taxon, species)
            if release_dir is None:
                return None
            ftp.cwd(release_dir)
            dirs = self._get_subdirs(ftp)

        # Manually filter suppressed directories
        if release_dir == "all_assembly_versions":
            dirs = [dir for dir in dirs if dir != "suppressed"]
        return dirs

    # ---- Genomic Asset Fetching ----
    @override
    def _matches_checksum(self, file_path: Path, expected_checksum: str) -> bool:
        """Checks whether a file maches the expected checksum from the md5checksums.txt file using an md5 hash.

        Arguments:
            file_path {pathlib.Path} -- Path to a local file.
            expected_checksum {str} -- Expected checksum of the file.

        Notes:
            The checksum type used is an md5 hash.

        Returns:
            bool -- Whether the file matches the expected checksum.
        """
        with open(file_path, "rb") as f:
            digest = hashlib.file_digest(f, "md5")

        return digest.hexdigest() == expected_checksum

    def fetch_genomic_entity(self, entity: NcbiGenomicEntity) -> dict[str, str]:
        """Fetches genomic entity from cache or downloads it if not cached yet.

        Arguments:
            entity {GenomicEntity} -- Genomic entity to fetch.

        Raises:
            RuntimeError: No cache_dir was set.
            RuntimeError: A resource cannot be located.
            RuntimeError: A download fails.
            RuntimeError: A checksum cannot be verified.

        Returns:
            dict[str, str] -- dict containing the file paths and the resolved release and assembly.
        """
        self.mode = self._get_mode(entity)
        self.mode._validate_taxon_and_source()

        fetched_genomic_entity = super().fetch_genomic_entity(entity)

        annotation_file = fetched_genomic_entity["annotation_file"]
        sequence_file = fetched_genomic_entity["sequence_file"]

        sequence_remote_dir, _, sequence_mapping_file = self._get_params("fasta")
        annotation_remote_dir, _, annotation_mapping_file = self._get_params("gtf")

        self._apply_mapping_chr_names("fasta", sequence_file, sequence_remote_dir, sequence_mapping_file)
        self._apply_mapping_chr_names("gtf", annotation_file, annotation_remote_dir, annotation_mapping_file)

        return {
            "annotation_file": annotation_file,
            "sequence_file": sequence_file,
            "annotation_release": fetched_genomic_entity["annotation_release"],
            "genome_assembly": fetched_genomic_entity["genome_assembly"],
        }

    @override
    def _get_entity_context(self, entity: NcbiGenomicEntity) -> GenomicEntityContext:
        """Retrieves context for a genomic entity.

        Arguments:
            entity {GenomicEntity} -- Genomic entity to get context for.

        Returns:
            GenomicEntityContext -- Context for the requested entity.
        """
        sequence_remote_dir, sequence_remote_filename, _ = self._get_params("fasta")
        annotation_remote_dir, annotation_remote_filename, _ = self._get_params("gtf")

        annotation_release = (
            self.annotation_release_from_ncbi
            if self.annotation_release_from_ncbi is not None
            else f"{self.assembly_accession}_{self.assembly_name}"
        )

        return GenomicEntityContext(
            annotation_remote_dir,
            annotation_remote_filename,
            sequence_remote_dir,
            sequence_remote_filename,
            annotation_release,
            self.assembly_name,
            self.assembly_accession,
        )

    def __init__(
        self,
        dir_output: str,
    ) -> None:
        """Constructor for the FtpLoaderNCBI class."""
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

    def _get_mode(self, entity: NcbiGenomicEntity) -> ModeStrategy:
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
                    f"mode '{self.mode}' is not supported. Supported values for source='ncbi' are: {allowed_modes}."
                )

    def _apply_mapping_chr_names(
        self,
        file_type: _TYPES_FILE,
        ftp_file: str,
        ftp_directory: str,
        ftp_mapping_file: str,
    ) -> None:
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
        All species annotated with the NCBI Eukaryotic Genome Annotation Pipeline should have a
        README_{annotation name} file that contains the annotation name used by NCBI. In case
        this is not found, keep the default value ("no_annotation_info").
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
