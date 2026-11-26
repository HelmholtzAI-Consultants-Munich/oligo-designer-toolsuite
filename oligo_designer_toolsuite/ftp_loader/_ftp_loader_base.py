############################################
# imports
############################################

import datetime
import gzip
import hashlib
import os
import re
import shutil
from abc import ABC, abstractmethod
from dataclasses import dataclass
from email.utils import formatdate, parsedate_to_datetime
from ftplib import FTP
from pathlib import Path
from typing import ClassVar, Generic, TypedDict, TypeVar, get_args

import requests

from oligo_designer_toolsuite._constants import _TYPES_FILE, _TYPES_FILE_SEQ

############################################
# FTP Classes
############################################


@dataclass(frozen=True)
class GenomicEntity:
    """Genomic entity for use with genomic databases."""

    taxon: str | None
    species: str
    release: str


@dataclass(frozen=True)
class RemoteInfo:
    dir: str
    name: str


@dataclass(frozen=True)
class GenomicEntityContext:
    """Context about a genomic entity for genomic databases."""

    file_remote_info: dict[_TYPES_FILE, RemoteInfo]
    annotation_release: str
    genome_assembly: str
    accession: str | None  # Ensembl doesn't use GCF/GCA accessions in filenames


class FetchedFiles(TypedDict):
    gtf: str | None
    gff: str | None
    fasta: str | None


class FetchedGenomicEntity(TypedDict):
    annotation_release: str
    genome_assembly: str
    fetched_files: dict[_TYPES_FILE, str]


class FtpUtil:
    def __init__(
        self,
        host: str,
        cache_dir: Path | None,
        name: str,
        download_chunk_size: int = 10 * 1024 * 1024,
    ):
        self.host = host
        self.cache_dir = cache_dir
        self.name = name
        self.DOWNLOAD_CHUNK_SIZE = download_chunk_size

    def get_exact_filename(self, ftp_directory: str, pattern: str) -> str:
        """
        Get exact filename by matching pattern against entries in FTP directory

        :param ftp_directory: The directory on the server where to look for the pattern
        :type ftp_directory: str
        :param pattern: The pattern that should be matched
        :type pattern: str

        """
        ftp = FTP(self.host)
        ftp.login()  # login to ftp server
        ftp.cwd(ftp_directory)  # move to directory

        files = ftp.nlst()
        exact_file_name = None

        for file in files:
            if re.match(pattern, file):
                exact_file_name = file

        ftp.quit()

        if exact_file_name is None:
            raise FileNotFoundError(f"File '{pattern}' not found on FTP server.")

        return exact_file_name

    def download(self, dir: str, remote_filename: str) -> Path:
        """
        Downloads a remote file if it changed.

        This function avoids redownloads if the file is already present and up
        to date. Since compressed files are deleted after extraction, this acts
        as a caching solution for genomic downloads.

        :param dir: Remote directory containing the file.
        :type dir: str
        :param remote_filename: Name of the remote file.
        :type remote_filename: str
        :return: Local file path of the downloaded resource.
        :rtype: pathlib.Path
        :raises RuntimeError: If no caching directory is set.
        """
        url = f"https://{self.host}/{dir}/{remote_filename}"
        url_hash = hashlib.md5(url.encode()).hexdigest()

        # Build local file path to save file at
        if self.cache_dir is None:
            raise RuntimeError("No caching directory set for genomic downloads.")
        file_path = (self.cache_dir / self.name / f"{url_hash}-{remote_filename}").resolve()
        file_path.parent.mkdir(parents=True, exist_ok=True)

        headers: dict[str, str] = {}

        # Check whether the file was downloaded before
        if file_path.exists():
            mtime = file_path.stat().st_mtime
            headers["if-modified-since"] = formatdate(mtime, usegmt=True)

        with requests.get(url, headers=headers, stream=True) as response:
            response.raise_for_status()

            if response.status_code == requests.codes.not_modified:
                # Local file is already up-to-date
                return file_path

            if response.status_code == requests.codes.ok:
                # Download the remote file
                with open(file_path, "wb") as f:
                    for chunk in response.iter_content(chunk_size=self.DOWNLOAD_CHUNK_SIZE):
                        f.write(chunk)

                # Set local modification time to remote's 'last-modified' to use for future requests
                if last_modified := response.headers.get("last-modified"):
                    new_mtime = parsedate_to_datetime(last_modified).timestamp()
                    os.utime(
                        file_path,
                        times=(datetime.datetime.now().timestamp(), new_mtime),
                    )
        return file_path


TGenomicEntity = TypeVar("TGenomicEntity")


class GenomicDataBaseConfigBase:
    name: ClassVar[str]
    host: ClassVar[str]
    base_path: ClassVar[str]
    checksums_filename: ClassVar[str]
    allowlist: set[str] | None


class GenomicDatabaseDownloaderBase(ABC, Generic[TGenomicEntity]):
    """
    Base class for downloading genomic data and postprocessing downloaded files.

    This class serves as the base interface for genomic databases such as Ensembl and NCBI and should not be used directly.

    For details on the caching procedure see "Caching FASTA Files" in the
    developer documentation.

    :param cache_dir: Path to the local directory used as a download cache.
    :type cache_dir: pathlib.Path | None
    :cvar str name: Name of the genomic database.
    :cvar str host: Hostname of the genomic database service.
    :cvar str base_path: Path to the subdirectory containing genomic files.
    :cvar str checksums_filename: Filename of checksum files for this database.
    """

    name: ClassVar[str]
    host: ClassVar[str]
    base_path: ClassVar[str]
    checksums_filename: ClassVar[str]

    def __init__(
        self,
        cache_dir: Path | None = None,
    ) -> None:
        """Constructor for the GenomicDatabaseDownloaderBase class."""
        self.cache_dir = cache_dir
        self.ftp_util = FtpUtil(self.host, self.cache_dir, self.name)

    @abstractmethod
    def _matches_checksum(self, file_path: Path, expected_checksum: str) -> bool:
        """
        Check whether a file matches the expected checksum.

        :param file_path: Path to a local file.
        :type file_path: pathlib.Path
        :param expected_checksum: Expected checksum of the file.
        :type expected_checksum: str
        :return: True if the file matches the checksum, otherwise False.
        :rtype: bool
        :raises NotImplementedError: If not implemented by subclasses.
        """

    def _download_and_process(self, dir: str, remote_filename: str, expected_checksum: str | None) -> Path:
        """
        Download a remote file and process it as needed.

        Processing includes optional checksum verification and extraction of gzip-compressed files.

        Since ``expected_checksum`` is part of the caching key and checksum files are not cached, a redownload will occur whenever the checksum changes.

        :param dir: Remote directory containing the to-be-downloaded file.
        :type dir: str
        :param remote_filename: Name of the remote file to download.
        :type remote_filename: str
        :param expected_checksum: Expected checksum of the file.
        :type expected_checksum: str | None
        :return: Local file path of the processed resource.
        :rtype: pathlib.Path
        :raises RuntimeError: If checksum verification fails
        """
        file_path = self.ftp_util.download(dir, remote_filename)

        # Verify checksum if provided. On mismatch, delete the cached file and retry the
        # download once before giving up.
        if expected_checksum is not None and not self._matches_checksum(file_path, expected_checksum):
            file_path.unlink(missing_ok=True)
            file_path = self.ftp_util.download(dir, remote_filename)
            if not self._matches_checksum(file_path, expected_checksum):
                raise RuntimeError(f"Checksum for {remote_filename} does not match.")

        # Extract file if extension is .gz
        if file_path.suffix == ".gz":
            file_path = self._decompress_gzip(file_path)

        return file_path

    @abstractmethod
    def _get_entity_context(
        self,
        entity: TGenomicEntity,
        files: list[_TYPES_FILE],
        sequence_nature: _TYPES_FILE_SEQ,
    ) -> GenomicEntityContext:
        """
        Retrieve context information for a genomic entity.

        :param entity: Genomic entity to get context for.
        :type entity: GenomicEntity
        :param files: List of file types to fetch.
        :type files: list[_TYPES_FILE]
        :param sequence_nature: Sequence Types to fetch.
        :type sequence_nature: _TYPES_FILE_SEQ
        :return: Context for the requested entity.
        :rtype: GenomicEntityContext
        """

    def _parse_checksum_line(self, line: str) -> tuple[str, str]:
        """
        Parse a line from a checksum file.

        Expects checksum lines to start with a checksum value and end with a
        filename, optionally containing additional whitespace-separated columns.
        If filenames start with ``./`` this prefix is removed.

        Subclasses may override this method if their checksum file format differs.

        :param line: Line from a checksum file.
        :type line: str
        :return: Tuple of ``(filename, checksum)``.
        :rtype: tuple[str, str]
        :raises ValueError: If the line cannot be parsed.
        """
        checksum, *_, filename = line.strip().split()
        filename = filename.removeprefix("./")
        return filename, checksum

    def _get_checksum_map(self, context: GenomicEntityContext) -> dict[str, str]:
        """
        Build a mapping of filenames to checksums.

        :param context: Context of the requested genomic entity.
        :type context: GenomicEntityContext
        :return: Mapping from filename to checksum.
        :rtype: dict[str, str]
        """
        # Combined checksum map for annotation (GTF) and sequence (FASTA) files
        filename_to_checksum_map: dict[str, str] = {}

        # Use a set for both directories -> if they are equal, the checksums are only downloaded once
        for dir in {remote_info.dir for remote_info in context.file_remote_info.values()}:
            # Always download without caching in case checksums changed
            checksums_path = self.ftp_util.download(dir, self.checksums_filename)
            with open(checksums_path) as checksums_file:
                for line in checksums_file:
                    try:
                        filename, checksum = self._parse_checksum_line(line)
                        filename_to_checksum_map[filename] = checksum
                    except ValueError:
                        pass

        return filename_to_checksum_map

    def _get_checksum(self, checksum_map: dict[str, str], filename: str) -> str:
        """
        Retrieve a checksum from a checksum map.

        Converts a ``KeyError`` into a ``RuntimeError`` for easier downstream
        handling.

        :param checksum_map: Mapping of filenames to checksums.
        :type checksum_map: dict[str, str]
        :param filename: Filename to look up.
        :type filename: str
        :return: Checksum of the specified file.
        :rtype: str
        :raises RuntimeError: If the file is not present in the checksum map.
        """
        try:
            return checksum_map[filename]
        except KeyError as e:
            raise RuntimeError(
                f"Required file missing in {self.name}'s {self.checksums_filename}: {e}"
            ) from e

    def fetch_genomic_entity(
        self,
        entity: TGenomicEntity,
        files: list[_TYPES_FILE] | None = None,
        sequence_nature: _TYPES_FILE_SEQ = "dna",
    ) -> FetchedGenomicEntity:
        """
        Fetch a genomic entity from cache or download it.

        :param entity: Genomic entity to fetch.
        :type entity: TGenomicEntity
        :param files: Requested file types.
        :type files: list[_TYPES_FILE] | None
        :param sequence_nature: Type of sequence data to retrieve.
        :type sequence_nature: _TYPES_FILE_SEQ
        :return: Dictionary containing file paths together with the resolved annotation release and genome assembly.
        :rtype: FetchedGenomicEntity
        :raises RuntimeError: If no cache directory is configured.
        :raises RuntimeError: If a resource cannot be located.
        :raises RuntimeError: If a download fails.
        :raises RuntimeError: If checksum verification fails.
        """

        if files is None:
            files = ["fasta", "gtf"]

        # Obtain GenomicEntityContext (file names, location, etc.) using subclass-specific implementation
        context = self._get_entity_context(entity, files, sequence_nature)
        # Obtain filename-to-checksum map
        checksum_map = self._get_checksum_map(context)

        fetched_files: dict[_TYPES_FILE, str] = {}

        for file_type, remote_info in context.file_remote_info.items():
            # Verify checksum is known
            checksum = self._get_checksum(checksum_map, remote_info.name)
            # Download from source, verify checksum of compressed file and extract if necessary
            annotation_file = self._download_and_process(remote_info.dir, remote_info.name, checksum)
            fetched_files[file_type] = str(annotation_file)

        return {
            "fetched_files": {**fetched_files},
            "annotation_release": context.annotation_release,
            "genome_assembly": context.genome_assembly,
        }

    def _decompress_gzip(self, file_path: Path) -> Path:
        """
        Decompress a gzip file and remove the original archive.

        :param file_path: Path to the gzip-compressed file.
        :type file_path: pathlib.Path
        :return: Path to the decompressed file.
        :rtype: pathlib.Path
        """

        extracted_file_path = file_path.with_suffix("")  # remove .gz extension
        with gzip.open(file_path, "rb") as archive:
            with open(extracted_file_path, "wb") as extract:
                shutil.copyfileobj(archive, extract)

            # Delete compressed file
            file_path.unlink()
            file_path = extracted_file_path

        return file_path

    def _check_file_type(self, file_type: _TYPES_FILE) -> None:
        """
        Check whether the provided file type is supported.

        :param file_type: File type to validate.
        :type file_type: _TYPES_FILE
        """
        options = get_args(_TYPES_FILE)
        assert file_type in options, f"File type not supported! '{file_type}' is not in {options}."

    def _check_sequence_nature_type(self, sequence_nature: _TYPES_FILE_SEQ) -> None:
        """
        Check whether the provided sequence nature is supported.

        :param sequence_nature: Sequence nature to validate.
        :type sequence_nature: _TYPES_FILE_SEQ
        """
        options = get_args(_TYPES_FILE_SEQ)
        assert (
            sequence_nature in options
        ), f"Sequence nature type not supported! '{sequence_nature}' is not in {options}."
