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
from typing import ClassVar, get_args

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
class GenomicEntityContext:
    """Context about a genomic entity for genomic databases."""

    annotation_remote_dir: str
    annotation_remote_filename: str
    sequence_remote_dir: str
    sequence_remote_filename: str
    annotation_release: str
    genome_assembly: str
    accession: str | None  # Ensembl doesn't use GCF/GCA accessions in filenames


class FtpLoaderBase(ABC):
    """
    A base class for downloading files via FTP and postprocessing the downloaded files.

    :param dir_output: Directory path where output files will be saved.
    :type dir_output: str
    """

    """The base class of genomic databases, which serve as our interfaces with public services like Ensembl and NCBI.

    Notes:
        This is a base class and should not be used directly.
        For details on the caching procedure see 'Caching FASTA Files' in the developer documentation.

    Class Variables:
        name {str} -- Name of the genomic database.
        host {str} -- Hostname of the genomic database service.
        base_path {str} -- Path to the subdirectory containing genomic files.
        checksums_filename {str} -- Filename of checksum files for this database.
    """

    name: ClassVar[str]
    host: ClassVar[str]
    base_path: ClassVar[str]
    checksums_filename: ClassVar[str]

    def __init__(
        self,
        cache_dir: Path | None = None,
        allowlist: set[str] | None = None,
        download_chunk_size: int = 10 * 1024 * 1024,
    ) -> None:
        """Initializes the genomic database.

        Keyword Arguments:
            cache_dir {pathlib.Path | None} -- Path to the local directory to use as a download cache. (default: {None})
            allowlist {set[str] | None} -- Set of subdirectories. If provided, queries will be limited to these. (default: {None})
        """
        self.cache_dir = cache_dir
        self.allowlist = allowlist
        self.DOWNLOAD_CHUNK_SIZE = download_chunk_size

    # ---- Directory Discovery ----
    def _get_subdirs(self, ftp: FTP) -> list[str]:
        """Retrieves a list of subdirectories at the current FTP cursor.

        Notes:
            This does not use ftplib.FTP.mlsd as Ensembl's FTP server does not support that option.

        Arguments:
            ftp {ftplib.FTP} -- Active FTP handler.

        Returns:
            list[str] -- list of subdirectories at the current FTP cursor.
        """
        file_lines: list[str] = []
        ftp.dir(file_lines.append)

        def parse_file_line(line: str) -> tuple[str, str]:
            """Parses filename and file permissions from a line returned by ftplib.FTP.dir.

            Arguments:
                line {str} -- Line returned from ftplib.FTP.dir.

            Returns:
                tuple[str, str] -- (filename, file permissions).
            """
            perms, *_, filename = line.split(maxsplit=8)
            # Normalize file permissions
            perms = perms.lower()
            # The filename might be a symbolic link (e.g. "current -> release-116")
            # Split this to just get the link's name (e.g. "current")
            filename = filename.split(maxsplit=1)[0]
            return filename, perms

        def is_likely_directory(filename: str, perms: str) -> bool:
            """Checks whether a file is likely a directory.

            Notes:
                Returns True for symbolic links if is does not have specific file extensions.
                As this is a very inaccurate heuristic, the function may misclassify symbolic links.

            Arguments:
                filename {str} -- Filename to check.
                perms {str} -- Permissions of the file.

            Returns:
                bool -- Whether the file is likely a directory.
            """
            non_dir_extensions = {"txt", "pdf", "gz"}
            _, _, file_suffix = filename.rpartition(".")

            return perms[0] == "d" or (perms[0] == "l" and file_suffix not in non_dir_extensions)

        dirs = [
            filename
            for filename, perms in map(parse_file_line, file_lines)
            if is_likely_directory(filename, perms)
        ]

        return sorted(dirs)  # sort for determinism

    def _filter_allowlist(self, dirs: list[str]) -> list[str]:
        """Filters the directories according to the allowlist.

        Arguments:
            dirs {list[str]} -- list of directories.

        Returns:
            list[str] -- list of all directories that are also in the allowlist.
        """
        if self.allowlist is not None:
            return sorted(list(set(dirs).intersection(self.allowlist)))  # sort for determinism
        return dirs

    def _get_subdirectories(self, dirs: list[str], ftp: FTP) -> dict[str, list[str]]:
        """Retrieves all subdirectories of all directories in dirs.

        Arguments:
            dirs {list[str]} -- list of directories.
            ftp {ftplib.FTP} -- Active FTP handler.

        Returns:
            dict[str, list[str]] -- dict mapping directories from dirs to their respective subdirectories.
        """
        subdirectories_by_directory: dict[str, list[str]] = {}
        for dir in dirs:
            _ = ftp.cwd(f"{self.base_path}/{dir}")
            subdirectories_by_directory[dir] = self._get_subdirs(ftp)
        return subdirectories_by_directory

    @abstractmethod
    def _build_species_mapping(self, species_dirs: dict[str, list[str]]) -> dict[str, list[str]]:
        """Builds the result of fetch_species_mapping.

        Arguments:
            species_dirs {dict[str, list[str]]} -- dict mapping directories their respective subdirectories.

        Notes:
            This is an abstract method and must be implemented by subclasses.

        Returns:
            dict[str, list[str]] -- dict mapping genomic categories to specific variants.
        """

    def fetch_species_mapping(self) -> dict[str, list[str]]:
        """Fetches all available species for this genomic database.

        Notes:
            For Ensembl, this returns a dict mapping species to their available annotation releases.
            For NCBI, this returns a dict mapping taxons to their respective species.
            This may be confusing and is unfortunate from a maintainability standpoint.

        Returns:
            dict[str, list[str]] -- dict mapping genomic categories to specific variants.
        """
        with FTP(self.host) as ftp:
            ftp.login()
            ftp.cwd(self.base_path)
            top_dirs = self._get_subdirs(ftp)
            top_dirs = self._filter_allowlist(top_dirs)
            species_dirs = self._get_subdirectories(top_dirs, ftp)
        return self._build_species_mapping(species_dirs)

    # ---- Genomic Asset Fetching ----
    def _download(self, dir: str, remote_filename: str) -> Path:
        """Downloads a remote file if it changed.

        Arguments:
            dir {str} -- Remote directory where the to-be-dowloaded file is located.
            remote_filename {str} -- Name of the remote file to download.

        Notes:
            This function avoids redownloads if the file is already present and up-to-date.
            Since compressed files are deleted after extraction, this does not suffice as
            a caching solution for genomic downloads - see _download_and_process for that.

        Raises:
            RuntimeError: No caching directory set for genomic downloads.

        Returns:
            pathlib.Path -- Local file path of the downloaded resource.
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
                    os.utime(file_path, times=(datetime.datetime.now().timestamp(), new_mtime))
        return file_path

    @abstractmethod
    def _matches_checksum(self, file_path: Path, expected_checksum: str) -> bool:
        """Checks whether a file maches the expected checksum.

        Arguments:
            file_path {pathlib.Path} -- Path to a local file.
            expected_checksum {str} -- Expected checksum of the file.

        Notes:
            This is an abstract method and must be implemented by subclasses.

        Returns:
            bool -- Whether the file matches the expected checksum.
        """

    def _download_and_process(self, dir: str, remote_filename: str, expected_checksum: str | None) -> Path:
        """Downloads a remote file and processes it as needed.

        Arguments:
            dir {str} -- Remote directory where the to-be-dowloaded file is located.
            remote_filename {str} -- Name of the remote file to download.
            expected_checksum {str | None} -- Expected checksum of the file.

        Notes:
            Processing includes optional verification of a checksum and extraction of .gz files.

            This function is decorated with our file cache to serve as the level 2 cache.
            Since expected_checksum is part of the caching key and the download of checksums
            is not being cached, a redownload will occur if the checksum changes.

        Raises:
            RuntimeError: Checksum for provided file does not match.

        Returns:
            pathlib.Path -- Local file path of the downloaded resource.
        """
        file_path = self._download(dir, remote_filename)

        # Verify checksum if provided. On mismatch, delete the cached file and retry the
        # download once before giving up.
        if expected_checksum is not None and not self._matches_checksum(file_path, expected_checksum):
            file_path.unlink(missing_ok=True)
            file_path = self._download(dir, remote_filename)
            if not self._matches_checksum(file_path, expected_checksum):
                raise RuntimeError(f"Checksum for {remote_filename} does not match.")

        # Extract file if extension is .gz
        if file_path.suffix == ".gz":
            file_path = self._decompress_gzip(file_path)

        return file_path

    @abstractmethod
    def _get_entity_context(self, entity: GenomicEntity) -> GenomicEntityContext:
        """Retrieves context for a genomic entity.

        Arguments:
            entity {GenomicEntity} -- Genomic entity to get context for.

        Notes:
            This is an abstract method and must be implemented by subclasses.

        Returns:
            GenomicEntityContext -- Context for the requested entity.
        """

    def _parse_checksum_line(self, line: str) -> tuple[str, str]:
        """Parses a line of a checksums file.

        Arguments:
            line {str} -- Line of a checksums file.

        Notes:
            Expects checksum lines to start with a checksum and end in a filename with
            an arbitrary amount of other data points inbetween, separated by whitespace.
            If filenames start with "./", this prefix will be removed.

            Subclasses may want to override this method for their specific checksums file format.

        Raises:
            ValueError: Parsing was unsuccessful.

        Returns:
            tuple[str, str] -- (filename, checksum).
        """
        checksum, *_, filename = line.strip().split()
        filename = filename.removeprefix("./")
        return filename, checksum

    def _get_checksum_map(self, context: GenomicEntityContext) -> dict[str, str]:
        """Returns map of filenames to their respective checksum.

        Arguments:
            context {GenomicEntityContext} -- Context for the requested entity.

        Returns:
            dict[str, str] -- dict mapping filenames to their respective checksum.
        """
        # Combined checksum map for annotation (GTF) and sequence (FASTA) files
        filename_to_checksum_map: dict[str, str] = {}

        # Use a set for both directories -> if they are equal, the checksums are only downloaded once
        for dir in {context.annotation_remote_dir, context.sequence_remote_dir}:
            # Always download without caching in case checksums changed
            checksums_path = self._download(dir, self.checksums_filename)
            with open(checksums_path) as checksums_file:
                for line in checksums_file:
                    try:
                        filename, checksum = self._parse_checksum_line(line)
                        filename_to_checksum_map[filename] = checksum
                    except ValueError:
                        pass

        return filename_to_checksum_map

    def _get_checksum(self, checksum_map: dict[str, str], filename: str) -> str:
        """Gets a checksum from a checksum map and handles key errors.

        Arguments:
            checksum_map {dict[str, str]} -- dict mapping filenames to their respective checksum.
            filename {str} -- Name of the file to look up in the checksum map.

        Notes:
            This method's purpose is to map the usual KeyError to a custom RuntimeError
            for easier downstream error catching and handling.

        Raises:
            RuntimeError: File not found in checksum map.

        Returns:
            str -- Checksum of the specified file.
        """
        try:
            return checksum_map[filename]
        except KeyError as e:
            raise RuntimeError(
                f"Required file missing in {self.name}'s {self.checksums_filename}: {e}"
            ) from e

    def fetch_genomic_entity(self, entity: GenomicEntity) -> dict[str, str]:
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

        # Obtain GenomicEntityContext (file names, location, etc.) using subclass-specific implementation
        context = self._get_entity_context(entity)
        # Obtain filename-to-checksum map
        checksum_map = self._get_checksum_map(context)

        # --- Annotation (GTF) ---
        # Verify checksum is known
        annotation_checksum = self._get_checksum(checksum_map, context.annotation_remote_filename)
        # Download from source, verify checksum of compressed file and extract if necessary
        annotation_file = self._download_and_process(
            context.annotation_remote_dir, context.annotation_remote_filename, annotation_checksum
        )

        # --- Sequence (FASTA) ---
        # Verify checksum is known
        sequence_checksum = self._get_checksum(checksum_map, context.sequence_remote_filename)
        # Download from source, verify checksum of compressed file and extract if necessary
        sequence_file = self._download_and_process(
            context.sequence_remote_dir, context.sequence_remote_filename, sequence_checksum
        )

        return {
            "annotation_file": str(annotation_file),
            "sequence_file": str(sequence_file),
            "annotation_release": context.annotation_release,
            "genome_assembly": context.genome_assembly,
        }

    def _decompress_gzip(self, file_path: Path) -> Path:
        """
        Decompresses a .gz file and removes the original compressed file.

        :param file_gzip: The path to the .gz file to be decompressed.
        :type file_gzip: str
        :return: The path to the decompressed file.
        :rtype: str
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

    def _get_exact_filename(self, ftp_directory: str, ftp_filename: str) -> str:
        ftp = FTP(self.host)
        ftp.login()  # login to ftp server
        ftp.cwd(ftp_directory)  # move to directory

        files = ftp.nlst()
        exact_file_name = None

        for file in files:
            if re.match(ftp_filename, file):
                exact_file_name = file

        ftp.quit()

        if exact_file_name is None:
            raise FileNotFoundError(f"File '{ftp_filename}' not found on FTP server.")

        return exact_file_name
