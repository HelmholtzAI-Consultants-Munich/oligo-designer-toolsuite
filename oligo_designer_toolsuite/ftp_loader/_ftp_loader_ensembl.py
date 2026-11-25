import os
import re
import subprocess
from pathlib import Path
from typing import ClassVar

from typing_extensions import override

from oligo_designer_toolsuite._constants import _TYPES_FILE, _TYPES_FILE_SEQ
from oligo_designer_toolsuite.ftp_loader._ftp_loader_base import (
    GenomicDataBaseConfigBase,
    GenomicDatabaseDownloaderBase,
    GenomicEntity,
    GenomicEntityContext,
    RemoteInfo,
)


class GenomicDataBaseConfigEnsembl(GenomicDataBaseConfigBase):
    name: ClassVar[str] = "ensembl"
    host: ClassVar[str] = "ftp.ensembl.org"
    base_path: ClassVar[str] = "/pub"
    checksums_filename: ClassVar[str] = "CHECKSUMS"


class GenomicDatabaseDownloaderEnsembl(GenomicDatabaseDownloaderBase[GenomicEntity]):
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

    name: ClassVar[str] = GenomicDataBaseConfigEnsembl.name
    host: ClassVar[str] = GenomicDataBaseConfigEnsembl.host
    base_path: ClassVar[str] = GenomicDataBaseConfigEnsembl.base_path
    checksums_filename: ClassVar[str] = GenomicDataBaseConfigEnsembl.checksums_filename

    def __init__(self, dir_output: str) -> None:
        """Constructor for the FtpLoaderEnsembl class."""
        super().__init__(Path(dir_output))
        self.assembly_name_placeholder = "[^\\.]*"

        self.file_type_folder = {"gff": "gff3", "gtf": "gtf", "fasta": "fasta"}

        self.file_type_ending = {
            "gff": "gff3.gz",
            "gtf": "gtf.gz",
            "fasta": "dna_sm.primary_assembly.fa.gz",  # soft-masked version of the genome
        }

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
    def _get_entity_context(
        self, entity: GenomicEntity, files: list[_TYPES_FILE], sequence_nature: _TYPES_FILE_SEQ
    ) -> GenomicEntityContext:
        """Retrieves context for a genomic entity.

        Arguments:
            entity {GenomicEntity} -- Genomic entity to get context for.

        Returns:
            GenomicEntityContext -- Context for the requested entity.
        """

        # resolve current to actual release
        annotation_release = self._get_annotation_release(entity.release)

        # Resolve filenames, remote directories and assembly
        remote_file_infos: dict[_TYPES_FILE, RemoteInfo] = {}
        for file_type in files:
            file_type_remote_dir, file_type_remote_filename = self._get_params(
                file_type, sequence_nature, entity.species, entity.release
            )
            remote_file_infos[file_type] = RemoteInfo(file_type_remote_dir, file_type_remote_filename)

        genome_assembly = self._get_genome_assembly(remote_file_infos)

        return GenomicEntityContext(
            remote_file_infos,
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

    def _get_params(
        self, file_type: _TYPES_FILE, sequence_nature: _TYPES_FILE_SEQ, species: str, annotation_release: str
    ) -> tuple[str, str]:
        if file_type == "fasta":
            return self._get_fasta_params(sequence_nature, species, annotation_release)
        else:
            return self._get_gtf_or_gff_params(species, annotation_release, file_type)

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

    def _get_gtf_or_gff_params(
        self, species: str, annotation_release: str, file_type: _TYPES_FILE
    ) -> tuple[str, str]:
        ftp_directory = f"pub/release-{annotation_release}/{self.file_type_folder[file_type]}/{species}/"
        ftp_file = f"{species.capitalize()}.{self.assembly_name_placeholder}.{annotation_release}.{self.file_type_ending[file_type]}"

        return ftp_directory, self.ftp_util.get_exact_filename(ftp_directory, ftp_file)

    def _parse_assembly_from_filename(self, file_name: str) -> str | None:
        """
        Try to parse assembly from annotation or sequence filename.

        Annotation file name examples:
            - Homo_sapiens.GRCh38.116.gtf.gz -> GRCh38
            - Acanthochromis_polyacanthus.ASM210954v1.116.gtf.gz -> ASM210954v1

        Sequence file name examples:
            - Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz -> GRCh38
            - Acanthochromis_polyacanthus.ASM210954v1.dna.toplevel.fa.gz -> ASM210954v1
        """

        match = re.search("\\.([^\\.]*)\\.", file_name)
        return match.group().replace(".", "") if match is not None else None

    def _get_genome_assembly(self, remote_file_infos: dict[_TYPES_FILE, RemoteInfo]) -> str:

        # Prefer assembly parsed from sequence; assembly parsed from annotation as fallback
        file_type_priority: list[_TYPES_FILE] = ["fasta", "gtf", "gff"]

        for file_type in file_type_priority:
            remote_file_info = remote_file_infos.get(file_type)
            if remote_file_info is None:
                continue
            assembly_from_annotation = self._parse_assembly_from_filename(remote_file_infos[file_type].name)
            if assembly_from_annotation is not None:
                return assembly_from_annotation

        return "unknown"
