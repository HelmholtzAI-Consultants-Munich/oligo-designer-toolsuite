import os
import re
import subprocess
from collections import defaultdict
from collections.abc import Callable
from ftplib import FTP
from pathlib import Path
from typing import ClassVar

from typing_extensions import override

from oligo_designer_toolsuite._constants import _TYPES_FILE_SEQ
from oligo_designer_toolsuite.ftp_loader._ftp_loader_base import (
    FtpLoaderBase,
    GenomicEntity,
    GenomicEntityContext,
)


class FtpLoaderEnsembl(FtpLoaderBase):
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

    """A genomic database to interface with Ensembl."""

    name: ClassVar[str] = "ensembl"
    host: ClassVar[str] = "ftp.ensembl.org"
    base_path: ClassVar[str] = "/pub"
    checksums_filename: ClassVar[str] = "CHECKSUMS"

    def __init__(self, dir_output: str, species: str, annotation_release: str) -> None:
        """Constructor for the FtpLoaderEnsembl class."""
        super().__init__(Path(dir_output))
        self.assembly_name_placeholder = "[^\\.]*"

        self.file_type_folder = {"gff": "gff3", "gtf": "gtf", "fasta": "fasta"}

        self.file_type_ending = {
            "gff": "gff3.gz",
            "gtf": "gtf.gz",
            "fasta": "dna_sm.primary_assembly.fa.gz",  # soft-masked version of the genome
        }

    # ---- Directory Discovery ----
    @override
    def _get_subdirectories(self, dirs: list[str], ftp: FTP) -> dict[str, list[str]]:
        """Retrieves all species for each release directory in dirs.

        Arguments:
            dirs {list[str]} -- list of release directories.
            ftp {ftplib.FTP} -- Active FTP handler.

        Returns:
            dict[str, list[str]] -- dict mapping releases from dirs to their respective species.
        """
        lookup_dirs = [f"{dir}/fasta" for dir in dirs]
        return super()._get_subdirectories(lookup_dirs, ftp)

    @override
    def _build_species_mapping(self, species_dirs: dict[str, list[str]]) -> dict[str, list[str]]:
        """Builds the result of fetch_species_mapping.

        Arguments:
            species_dirs {list[str]} -- dict mapping release directories to their respective species.

        Returns:
            dict[str, list[str]] -- dict mapping species to their available annotation releases.
        """

        def format_release_name(release_name: str) -> str:
            """Formats the release naming such that numeric releases are reduced to the plain number.

            Arguments:
                release_name {str} -- Name of the release to format.

            Returns:
                str -- Formatted release name.
            """
            return release_name.removeprefix("release-").removesuffix("/fasta")

        def reverse_list_dict(
            original_dict: dict[str, list[str]], key_formatter: Callable[[str], str]
        ) -> dict[str, list[str]]:
            """Reverses a dict of lists such that each list element is mapped to a list of keys it's contained in.

            Arguments:
                original_dict {dict[str, list[str]]} -- dict of lists.
                key_formatter {Callable[[str], str]} -- Function that formats keys before they are added to a result list.

            Returns:
                dict[str, list[str]] -- Reversed dict.
            """
            reversed_dict = defaultdict(list)
            for key, value_list in original_dict.items():
                for value in value_list:
                    reversed_dict[value].append(key_formatter(key))
            return dict(reversed_dict)

        # Turn version->list-of-species map into species->list-of-versions map
        return reverse_list_dict(species_dirs, format_release_name)

    # ---- Genomic Asset Fetching ----
    @override
    def _matches_checksum(self, file_path: Path, expected_checksum: str) -> bool:
        """Checks whether a file maches the expected checksum from the CHECKSUMS file using the 'sum' command.

        Arguments:
            file_path {pathlib.Path} -- Path to a local file.
            expected_checksum {str} -- Expected checksum of the file.

        Notes:
            The checksum type used is a 16-bit BSD checksum. There does not seem to be a
            neater Python implementation for this than simply calling the `sum` command.

        Returns:
            bool -- Whether the file matches the expected checksum.
        """
        try:
            result = subprocess.run(["sum", file_path], capture_output=True, check=True, text=True)
            computed_checksum = result.stdout.split(maxsplit=1)[0]
            return computed_checksum == expected_checksum
        except subprocess.CalledProcessError:
            return False

    @override
    def _get_entity_context(self, entity: GenomicEntity) -> GenomicEntityContext:
        """Retrieves context for a genomic entity.

        Arguments:
            entity {GenomicEntity} -- Genomic entity to get context for.

        Returns:
            GenomicEntityContext -- Context for the requested entity.
        """

        # resolve current to actual release
        annotation_release = self._get_annotation_release(entity.release)

        # Resolve filenames, remote directories and assembly
        sequence_remote_dir, sequence_remote_filename = self._get_fasta_params(
            "dna", entity.species, annotation_release
        )
        annotation_remote_dir, annotation_remote_filename = self._get_gtf_or_gff_params(
            entity.species, annotation_release
        )

        genome_assembly = self._get_genome_assembly(annotation_remote_filename, sequence_remote_filename)

        return GenomicEntityContext(
            annotation_remote_dir,
            annotation_remote_filename,
            sequence_remote_dir,
            sequence_remote_filename,
            annotation_release,
            genome_assembly,
            accession=None,
        )

    def _get_annotation_release(self, annotation_release: str) -> str:
        """
        Constructs the FTP directory path and file name based on the file type and sequence nature.

        :param annotation_release: The user set annotation release, which is either a number or 'current'
        :type annotation_release: str
        :return: The annotation release that can be downloaded from Ensembl.
        :rtype: str
        """
        if annotation_release == "current":
            file_version = self.ftp_util.download(
                "pub",
                "VERSION",
            )
            with open(file_version, "r") as handle:
                annotation_release = handle.readline().strip()
            os.remove(file_version)

        return annotation_release

    def _get_fasta_params(
        self, sequence_nature: _TYPES_FILE_SEQ, species: str, annotation_release: str
    ) -> tuple[str, str]:
        ftp_directory = (
            f"pub/release-{annotation_release}/{self.file_type_folder['fasta']}/{species}/{sequence_nature}/"
        )
        if sequence_nature == "dna":
            ftp_file = (
                f"{species.capitalize()}.{self.assembly_name_placeholder}.{self.file_type_ending['fasta']}"
            )
        else:
            ftp_file = f"{species.capitalize()}.{self.assembly_name_placeholder}.{sequence_nature}.fa.gz"

        return ftp_directory, self.ftp_util.get_exact_filename(ftp_directory, ftp_file)

    def _get_gtf_or_gff_params(self, species: str, annotation_release: str) -> tuple[str, str]:
        # TODO: how to handle missing file types here? Passthrough?
        ftp_directory = f"pub/release-{annotation_release}/{self.file_type_folder['gtf']}/{species}/"
        ftp_file = f"{species.capitalize()}.{self.assembly_name_placeholder}.{annotation_release}.{self.file_type_ending['gtf']}"

        return ftp_directory, self.ftp_util.get_exact_filename(ftp_directory, ftp_file)

    def _parse_assembly_from_filename(self, file_name: str) -> str | None:
        match = re.search("\\.([^\\.]*)\\.", file_name)
        return match.group().replace(".", "") if match is not None else None

    def _get_genome_assembly(self, annotation_filename: str, sequence_filename: str) -> str:
        # Try to parse assembly from annotation filename
        # e.g. Homo_sapiens.GRCh38.116.gtf.gz -> GRCh38
        # e.g. Acanthochromis_polyacanthus.ASM210954v1.116.gtf.gz -> ASM210954v1
        assembly_from_annotation = self._parse_assembly_from_filename(annotation_filename)

        # Try to parse assembly from sequence filename
        # e.g. Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz -> GRCh38
        # e.g. Acanthochromis_polyacanthus.ASM210954v1.dna.toplevel.fa.gz -> ASM210954v1
        assembly_from_sequence = self._parse_assembly_from_filename(sequence_filename)

        # Prefer assembly parsed from sequence; assembly parsed from annotation as fallback
        return assembly_from_sequence or assembly_from_annotation or "unknown"
