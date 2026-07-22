"""
Parse annotation, sequence, variant, and interval files used across the toolsuite.

GFF/GTF annotations, FASTA sequences, VCF variants, and BED intervals are the main
file types for region extraction, probe design, and coordinate conversion.
See :class:`GffParser`, :class:`FastaParser`, :class:`VCFParser`, and :class:`BedParser`.
"""

############################################
# imports
############################################

import gzip
import os
import pickle
import re
import subprocess
from typing import Any

import pandas as pd
from Bio import SeqIO
from cyvcf2 import VCF, Writer

from oligo_designer_toolsuite._constants import (
    SEPARATOR_FASTA_HEADER_FIELDS,
    SEPARATOR_FASTA_HEADER_FIELDS_LIST,
)
from oligo_designer_toolsuite._exceptions import ConfigurationError, FileFormatError

from ._checkers_and_helpers import cast_to_list

############################################
# GFF Parser Class
############################################


class GffParser:
    """
    Parse genome annotations from GFF or GTF files.

    GFF/GTF files describe features such as genes, transcripts, and exons on a
    reference sequence. Feature starts are 1-based. Use this class to check format,
    load annotations into a table, and optionally cache them as a pickle for reuse.
    """

    def __init__(self) -> None:
        """Constructor for the GffParser class."""
        self.GFF_HEADER = [
            "seqid",
            "source",
            "type",
            "start",
            "end",
            "score",
            "strand",
            "phase",
        ]
        self.R_SEMICOLON = re.compile(r"\s*;\s*")
        self.R_COMMA = re.compile(r"\s*,\s*")
        self.R_KEYVALUE = re.compile(r"(\s+|\s*=\s*)")

    def check_gff_format(self, file: str) -> bool:
        """
        Check whether a GFF/GTF annotation file exists and has valid content.

        Reads a sample of feature rows to confirm the file can be parsed. Useful
        before loading a full annotation into a pipeline step.

        :param file: Path to the GFF or GTF file to check.
        :type file: str
        :raises FileFormatError: If the file is missing or not a valid GFF/GTF.
        :return: ``True`` when the file exists and is a valid GFF/GTF.
        :rtype: bool
        """

        def _check_gff_content(file: str) -> bool:
            """
            Check whether the annotation file yields at least one feature row.

            Parses a short sample of the file. An empty result means the content
            is not usable as GFF/GTF.

            :param file: Path to the GFF or GTF file to sample.
            :type file: str
            :return: ``True`` if at least one feature row is found.
            :rtype: bool
            """
            gtf = self.parse_annotation_from_gff(file, target_lines=100)
            return any(gtf)

        if os.path.exists(file):
            if not _check_gff_content(file):
                raise FileFormatError(f"GFF file '{file}' has incorrect format. Expected valid GFF format.")
            return True
        else:
            raise FileFormatError(f"GFF file '{file}' does not exist.")

    def parse_annotation_from_gff(
        self,
        annotation_file: str,
        file_pickle: str | None = None,
        chunk_size: int = 10000,
        target_lines: int = 100000000,
    ) -> str | pd.DataFrame:
        """
        Parse a GFF/GTF annotation into a table of genomic features.

        Combines the standard GFF columns with attribute fields into one DataFrame.
        Feature starts are 1-based. Optionally write the table to a pickle for faster
        reloads later.

        :param annotation_file: Path to the GFF or GTF annotation file.
        :type annotation_file: str
        :param file_pickle: Optional path to save the table as a pickle. If set, that
            path is returned instead of the DataFrame.
        :type file_pickle: str | None
        :param chunk_size: Number of lines to process per batch.
        :type chunk_size: int
        :param target_lines: Maximum number of lines to read (useful for sampling).
        :type target_lines: int
        :return: Feature table, or the pickle path when ``file_pickle`` is set.
        :rtype: str | pd.DataFrame
        """
        csv_file, extra_info_file = self._split_annotation(
            annotation_file, chunk_size=chunk_size, target_lines=target_lines
        )

        info_df = self._info_to_df(extra_info_file, chunk_size=chunk_size)
        csv_df = pd.read_csv(
            csv_file, sep="\t", names=self.GFF_HEADER, header=None, dtype={"seqid": "string"}
        )

        csv_df.reset_index(inplace=True, drop=True)
        info_df.reset_index(inplace=True, drop=True)

        dataframe = pd.concat([csv_df, info_df], axis=1)

        os.remove(csv_file)
        os.remove(extra_info_file)

        if file_pickle:
            with open(file_pickle, "wb") as handle:
                pickle.dump(dataframe, handle)
            return file_pickle
        else:
            return dataframe

    def load_annotation_from_pickle(self, file_pickel: str) -> pd.DataFrame:
        """
        Load a previously saved GFF/GTF annotation table from a pickle file.

        Use this after :meth:`parse_annotation_from_gff` wrote a pickle, to avoid
        re-parsing the original annotation.

        :param file_pickel: Path to the pickle file with the annotation table.
        :type file_pickel: str
        :return: Table of genomic features from the pickle.
        :rtype: pd.DataFrame
        """
        dataframe_gff: pd.DataFrame = pickle.load(open(file_pickel, "rb"))

        return dataframe_gff

    def _split_annotation(self, annotation_file: str, chunk_size: int, target_lines: int) -> tuple[str, str]:
        """
        Split a GFF/GTF file into standard columns and attribute fields.

        Writes two temporary files: one with the eight core GFF columns, and one
        with the ninth-column attributes. Comment lines starting with ``#`` are skipped.

        :param annotation_file: Path to the GFF or GTF annotation file.
        :type annotation_file: str
        :param chunk_size: Number of lines to read per batch.
        :type chunk_size: int
        :param target_lines: Maximum number of lines to process.
        :type target_lines: int
        :return: Paths to the core-column file and the attributes file.
        :rtype: tuple[str, str]
        """
        csv_file = ".".join(annotation_file.split(".")[:-1]) + ".csv"
        extra_info_file = ".".join(annotation_file.split(".")[:-1]) + ".txt"

        finished = False
        lines_read = 0

        fn_open: Any = gzip.open if annotation_file.endswith(".gz") else open
        with fn_open(annotation_file, "r") as input_file:
            with open(csv_file, "w") as out_csv:
                with open(extra_info_file, "w") as out_extra_info:
                    while not finished and lines_read < target_lines:
                        csv_content_chunck = ""
                        extra_info_content_chunck = ""
                        for _ in range(chunk_size):
                            lines_read += 1
                            if lines_read > target_lines:
                                break
                            try:
                                line = next(input_file)
                                if not line.startswith("#"):
                                    csv_content_chunck += "\t".join(line.split("\t")[:8]) + "\n"
                                    extra_info_content_chunck += "\t".join(line.split("\t")[8:])
                            except:
                                finished = True

                        out_csv.write(csv_content_chunck)
                        out_extra_info.write(extra_info_content_chunck)

        return csv_file, extra_info_file

    def _parse_fields(self, line: str) -> dict[str, Any]:
        """
        Parse attribute key-value pairs from a GFF/GTF attributes column.

        Turns semicolon-separated tags such as gene and transcript IDs into a
        dictionary for joining with the core feature columns.

        :param line: Attributes text from one GFF/GTF feature row.
        :type line: str
        :return: Attribute names mapped to their values.
        :rtype: dict[str, Any]
        """
        result = {}

        # INFO field consists of "key1=value;key2=value;...".
        infos = [x for x in re.split(self.R_SEMICOLON, line) if x.strip()]

        for i, info in enumerate(infos, 1):
            # It should be key="value".
            try:
                key, _, value = re.split(self.R_KEYVALUE, info, maxsplit=1)
                key = key.strip("\"'")
                if key in self.GFF_HEADER:
                    key = f"attribe_{key}"
            # But sometimes it is just "value".
            except ValueError:
                key = f"INFO{i}"
                value = info
            # Ignore the field if there is no value.
            if value:
                result[key] = self._get_value(value)

        return result

    def _get_value(self, value: str) -> str | list[str] | None:
        """
        Normalize one attribute value from a GFF/GTF attributes field.

        Strips quotes, splits comma-separated lists, and treats empty, ``.``, and
        ``NA`` as missing.

        :param value: Raw attribute value text.
        :type value: str
        :return: Cleaned string, list of strings, or ``None`` if missing.
        :rtype: str | list[str] | None
        """
        if not value:
            return None

        # Strip double and single quotes and new lines.
        value = value.strip("\"'\n")

        # Return a list if the value has a comma.
        if "," in value:
            return re.split(self.R_COMMA, value)
        # These values are equivalent to None.
        elif value in ["", ".", "NA"]:
            return None
        else:
            return value

    def _info_to_df_chunk(self, data_chunk: list[str]) -> pd.DataFrame:
        """
        Convert a batch of GFF/GTF attribute lines into a table.

        Each line becomes one row of parsed attribute columns for later merge
        with the core GFF columns.

        :param data_chunk: Attribute lines from the temporary info file.
        :type data_chunk: list[str]
        :return: Table of parsed attributes for the batch.
        :rtype: pd.DataFrame
        """
        parsed_fields = [self._parse_fields(line) for line in data_chunk]
        return pd.DataFrame(parsed_fields)

    def _info_to_df(self, info_file: str, chunk_size: int) -> pd.DataFrame:
        """
        Load all GFF/GTF attribute lines from a temporary file into one table.

        Reads the file in batches so large annotations stay manageable, then
        concatenates the batches.

        :param info_file: Path to the temporary attributes file.
        :type info_file: str
        :param chunk_size: Number of lines to parse per batch.
        :type chunk_size: int
        :return: Combined table of attribute columns.
        :rtype: pd.DataFrame
        """
        info_dfs = []
        with open(info_file, "r") as info_f:
            data = info_f.readlines()
            n_lines = len(data)
            for i in range(0, n_lines, chunk_size):
                data_chunk = data[i : i + chunk_size]
                info_dfs.append(self._info_to_df_chunk(data_chunk))
        return pd.concat(info_dfs)


############################################
# Fasta Parser Class
############################################


class FastaParser:
    """
    Read, write, and merge nucleic acid sequences in FASTA format.

    FASTA files store reference genomes, extracted regions, and oligo sequences
    used throughout the toolsuite. Headers from the region generator carry a region
    ID plus optional metadata and 1-based coordinates
    (``chromosome:start-end(strand)``).
    """

    def __init__(self) -> None:
        """Constructor for the FastaParser class."""

    def check_fasta_format(self, file: str) -> bool:
        """
        Check whether a FASTA file exists and contains valid sequence records.

        Confirms the file can be opened as FASTA before downstream read or merge
        steps rely on it.

        :param file: Path to the FASTA file to check.
        :type file: str
        :raises FileFormatError: If the file is missing or not a valid FASTA.
        :return: ``True`` when the file exists and is a valid FASTA.
        :rtype: bool
        """

        def _check_fasta_content(file: str) -> bool:
            """
            Check whether the FASTA file yields at least one sequence record.

            An empty index means the content is not usable as FASTA.

            :param file: Path to the FASTA file to sample.
            :type file: str
            :return: ``True`` if at least one sequence record is found.
            :rtype: bool
            """
            fasta = SeqIO.index(file, "fasta")
            return any(fasta)  # False when `fasta` is empty, i.e. wasn't a FASTA file

        if os.path.exists(file):
            if not _check_fasta_content(file):
                raise FileFormatError(
                    f"FASTA file '{file}' has incorrect format. Expected valid FASTA format."
                )
            else:
                return True
        else:
            raise FileFormatError(f"FASTA file '{file}' does not exist.")

    def is_coordinate(self, entry: str) -> bool:
        """
        Check whether a header field looks like a genomic coordinate string.

        Matches the region-generator form ``chromosome:start-end(strand)``, where
        starts are 1-based.

        :param entry: Header field text to test.
        :type entry: str
        :return: ``True`` if the text matches the coordinate pattern.
        :rtype: bool
        """
        pattern = r"(\S+):(\d+)-(\d+)\(.*\)"
        return bool(re.match(pattern, entry))

    def parse_number(self, s: str) -> int | float | str:
        """
        Parse a string as an integer or float when possible.

        Used when reading numeric metadata from FASTA header key-value fields.
        Non-numeric text is returned unchanged.

        Examples::

            >>> parse_number("42")
            42
            >>> parse_number("3.14")
            3.14
            >>> parse_number("abc")
            "abc"
            >>> parse_number("5.0")
            5.0

        :param s: Text to parse.
        :type s: str
        :return: Integer, float, or the original string.
        :rtype: int | float | str
        """
        try:
            return int(s)
        except ValueError:
            pass

        try:
            return float(s)
        except ValueError:
            pass

        return s

    def get_fasta_regions(self, file_fasta_in: str) -> list[str]:
        """
        List unique region IDs from a FASTA file.

        Region IDs are taken from the first field of each sequence header, as
        produced by the genomic region generator.

        :param file_fasta_in: Path to the input FASTA file.
        :type file_fasta_in: str
        :return: Unique region identifiers found in the file.
        :rtype: list[str]
        """
        region_ids = []
        # use index instead of parse function for memory efficiency
        for idx in SeqIO.index(file_fasta_in, "fasta"):
            region, _, _ = self.parse_fasta_header(idx, parse_additional_info=False)
            region_ids.append(region)

        return list(set(region_ids))

    def read_fasta_sequences(
        self, file_fasta_in: str, region_ids: str | list[str] | None = None
    ) -> list[Any]:
        """
        Load sequences from a FASTA file, optionally filtered by region ID.

        When ``region_ids`` is set, only sequences whose header region matches are
        returned. Otherwise every record in the file is loaded.

        :param file_fasta_in: Path to the input FASTA file.
        :type file_fasta_in: str
        :param region_ids: One region ID, a list of IDs, or ``None`` for all regions.
        :type region_ids: str | list[str] | None
        :return: Sequence records from the file (filtered when IDs are given).
        :rtype: list[Any]
        """
        region_ids_set = set(region_ids) if region_ids else None

        # Open the file once and parse sequences accordingly
        if region_ids_set:
            fasta_sequences = []
            # use index instead of parse function for memory efficiency
            seq_record = SeqIO.index(file_fasta_in, "fasta")
            for idx in seq_record:
                region, _, _ = self.parse_fasta_header(
                    idx, parse_coordinates=False, parse_additional_info=False
                )
                if region in region_ids_set:
                    fasta_sequences.append(seq_record[idx])

        else:
            fasta_sequences = list(SeqIO.parse(file_fasta_in, "fasta"))

        return fasta_sequences

    def parse_fasta_header(
        self, header: str, parse_coordinates: bool = True, parse_additional_info: bool = True
    ) -> tuple[str, dict[str, list[Any]] | str, dict[str, list[Any]]]:
        """
        Parse region ID, coordinates, and metadata from a FASTA header.

        Headers from the region generator look like
        ``region_id::key=value::chr:start-end(strand)``. Coordinate starts are
        1-based. Set the parse flags to skip coordinates or metadata when not needed.

        :param header: FASTA header text (with or without a leading ``>``).
        :type header: str
        :param parse_coordinates: If ``True``, extract chromosome, start, end, and strand.
        :type parse_coordinates: bool
        :param parse_additional_info: If ``True``, parse key-value metadata into a dict;
            if ``False``, return that field as a string when present.
        :type parse_additional_info: bool
        :return: Region ID, additional info (dict or string), and coordinates dict.
        :rtype: tuple[str, dict[str, list[Any]] | str, dict[str, list[Any]]]
        """
        region: str = ""
        additional_info: dict[str, list[Any]] | str = {}
        coordinates: dict[str, list[Any]] = {
            "chromosome": [None],
            "start": [None],
            "end": [None],
            "strand": [None],
        }

        for header_entry in header.split(SEPARATOR_FASTA_HEADER_FIELDS):
            header_entry = header_entry.strip()
            if region == "":
                region = header_entry

            elif self.is_coordinate(header_entry):
                if parse_coordinates:
                    header_coordinates = header_entry.split(SEPARATOR_FASTA_HEADER_FIELDS_LIST)
                    coordinates = {}
                    for header_coordinate in header_coordinates:
                        coordinates.setdefault("chromosome", []).append(header_coordinate.split(":")[0])
                        coordinates.setdefault("start", []).append(
                            int(header_coordinate.split(":")[1].split("-")[0])
                        )
                        coordinates.setdefault("end", []).append(
                            int(header_coordinate.split(":")[1].split("-")[1].split("(")[0])
                        )
                        coordinates.setdefault("strand", []).append(
                            header_coordinate.split("(")[1].split(")")[0]
                        )

            else:
                info = header_entry
                # the additional info field should be parsed, save information in dict
                if parse_additional_info:
                    additional_info = {}
                    if SEPARATOR_FASTA_HEADER_FIELDS_LIST in info:
                        info_list = info.split(SEPARATOR_FASTA_HEADER_FIELDS_LIST)
                        for infos in info_list:
                            key, value = infos.split("=")
                            value_parsed = self.parse_number(value)
                            if key in additional_info:
                                additional_info[key].append(value_parsed)
                            else:
                                additional_info[key] = [value_parsed]

                    else:
                        key, value = info.split("=")
                        value_parsed = self.parse_number(value)
                        additional_info[key] = [value_parsed]
                else:
                    additional_info = str(info)

        return region, additional_info, coordinates

    def write_fasta_sequences(self, fasta_sequences: list[Any], file_out: str) -> None:
        """
        Write sequence records to a FASTA file.

        Overwrites ``file_out`` with the given records in standard FASTA layout.

        :param fasta_sequences: Sequence records to write.
        :type fasta_sequences: list[Any]
        :param file_out: Path of the output FASTA file.
        :type file_out: str
        :return: None
        """
        with open(file_out, "w") as handle_fasta:
            SeqIO.write(fasta_sequences, handle_fasta, "fasta")

    def merge_fasta_files(self, files_in: list[str], file_out: str, overwrite: bool = False) -> None:
        """
        Merge several FASTA files into one output file.

        Concatenates all sequence records from the inputs. Set ``overwrite`` to
        replace an existing output; otherwise records are appended.

        :param files_in: Paths of the FASTA files to merge.
        :type files_in: list[str]
        :param file_out: Path of the merged FASTA file.
        :type file_out: str
        :param overwrite: If ``True``, replace ``file_out``; if ``False``, append.
        :type overwrite: bool
        :return: None
        """
        files_in = cast_to_list(files_in)
        file_out_mode = "w" if overwrite else "a"
        with open(file_out, file_out_mode) as handle_fasta:
            for file_in in files_in:
                seq_record = SeqIO.index(file_in, "fasta")
                for idx in seq_record:
                    SeqIO.write(seq_record[idx], handle_fasta, "fasta")

    def index_fasta_file(self, file_fasta: str) -> None:
        """
        Create or refresh a ``.fai`` index for a FASTA file.

        Removes any existing index, then builds a new one with ``samtools faidx``.
        Call this after rewriting a FASTA so tools such as bedtools see matching
        sequence offsets.

        :param file_fasta: Path to the FASTA file to index.
        :type file_fasta: str
        :raises FileFormatError: If the FASTA is missing or indexing fails.
        :return: None
        """
        if not os.path.exists(file_fasta):
            raise FileFormatError(f"FASTA file '{file_fasta}' does not exist.")

        index_file = f"{file_fasta}.fai"

        # Remove existing index
        if os.path.exists(index_file):
            os.remove(index_file)

        # Use samtools faidx to create the index (same tool that bedtools getfasta uses)
        args = ["samtools", "faidx", file_fasta]
        process = subprocess.run(args)

        if process.returncode != 0:
            raise FileFormatError(
                f"Failed to create FASTA index file '{index_file}' using samtools faidx. "
                f"Please ensure samtools is installed and available in your PATH."
            )

        # Verify the index file was created
        if not os.path.exists(index_file):
            raise FileFormatError(f"Failed to create FASTA index file '{index_file}'.")


############################################
# VCF Parser Class
############################################


class VCFParser:
    """
    Read, write, and merge variant records in VCF format.

    VCF files store SNPs and other variants relative to a reference genome.
    Use this class to check format, load variants for filtering or annotation,
    write selected records, and merge multiple VCF inputs into one file.
    """

    def __init__(self) -> None:
        """Constructor for the VCFParser class."""

    def check_vcf_format(self, file: str) -> bool:
        """
        Check whether a VCF file exists and contains valid variant records.

        Confirms the file can be opened as VCF before read or merge steps use it.

        :param file: Path to the VCF file to check.
        :type file: str
        :raises FileFormatError: If the file is missing or not a valid VCF.
        :return: ``True`` when the file exists and is a valid VCF.
        :rtype: bool
        """

        def _check_vcf_content(file: str) -> bool:
            """
            Check whether the VCF file yields at least one variant record.

            An empty result means the content is not usable as VCF.

            :param file: Path to the VCF file to sample.
            :type file: str
            :return: ``True`` if at least one variant record is found.
            :rtype: bool
            """
            vcf = VCF(file)
            return any(vcf)  # False when `vcf` is empty, i.e. wasn't a vcf file

        if os.path.exists(file):
            if not _check_vcf_content(file):
                raise FileFormatError(f"VCF file '{file}' has incorrect format. Expected valid VCF format.")
            else:
                return True
        else:
            raise FileFormatError(f"VCF file '{file}' does not exist.")

    def read_vcf_variants(self, file: str) -> tuple[list[Any], Any]:
        """
        Load variant records from a VCF file.

        Checks format first, then returns every variant plus the VCF handle used
        for writing compatible output later.

        :param file: Path to the VCF file to read.
        :type file: str
        :raises FileFormatError: If the file is missing or not a valid VCF.
        :return: List of variant records and the VCF file handle.
        :rtype: tuple[list[Any], Any]
        """
        self.check_vcf_format(file)

        variants = []

        vcf_in = VCF(file)
        variants.extend(list(vcf_in))
        vcf_in.close()

        return variants, vcf_in

    def write_vcf_variants(self, vcf_variants: list, vcf_in: str, file_out: str) -> None:
        """
        Write variant records to a VCF file.

        Uses the input VCF handle for header and formatting so the output stays
        compatible with the source file.

        :param vcf_variants: Variant records to write.
        :type vcf_variants: list
        :param vcf_in: VCF handle that supplies header and format information.
        :type vcf_in: str
        :param file_out: Path of the output VCF file.
        :type file_out: str
        :return: None
        """
        vcf_out = Writer(file_out, vcf_in)
        for variant in vcf_variants:
            vcf_out.write_record(variant)
        vcf_out.close()

    def merge_vcf_files(self, files_in: list[str], file_out: str) -> None:
        """
        Merge several VCF files into one output file.

        Compresses, sorts, and indexes inputs as needed, then merges them with
        bcftools into a single uncompressed VCF.

        :param files_in: Paths of the VCF files to merge.
        :type files_in: list[str]
        :param file_out: Path of the merged VCF file.
        :type file_out: str
        :return: None
        """
        # Track compressed files that need to be cleaned up
        compressed_files_to_cleanup: list[str] = []

        args = ["bcftools", "merge", "--force-single", "--output-type", "v", "-o", file_out]
        for file_vcf in files_in:
            _, ext = os.path.splitext(file_vcf)
            if ext == ".vcf":
                file_vcf_compressed = f"{file_vcf}.gz"
                compressed_files_to_cleanup.append(file_vcf_compressed)

                args_compress = ["bcftools", "view", "-O", "z", "-o", file_vcf_compressed, file_vcf]
                subprocess.run(args_compress, check=True, stdout=subprocess.DEVNULL)

                file_vcf = file_vcf_compressed

            args_sort = ["bcftools", "sort", file_vcf, "-Oz", "-o", file_vcf]
            subprocess.run(args_sort, check=True, stdout=subprocess.DEVNULL)

            args_index = ["bcftools", "index", "-f", file_vcf]
            subprocess.run(args_index, check=True, stdout=subprocess.DEVNULL)

            args.append(file_vcf)

        subprocess.run(args, check=True, stdout=subprocess.DEVNULL)

        # Clean up compressed files and their index files
        for file_vcf_compressed in compressed_files_to_cleanup:
            if os.path.exists(file_vcf_compressed):
                os.remove(file_vcf_compressed)
            # Remove index files for compressed files
            for index_ext in [".csi", ".tbi"]:
                index_file = f"{file_vcf_compressed}{index_ext}"
                if os.path.exists(index_file):
                    os.remove(index_file)


############################################
# BED Parser Class
############################################


class BedParser:
    """
    Read and write genomic intervals in BED format.

    BED starts are 0-based. Region-generator FASTA headers and GFF/GTF annotations
    use 1-based starts; convert with :meth:`convert_start` when moving between those
    coordinate systems. On read, UCSC ``track``, ``browser``, and ``#`` lines are skipped.
    """

    BED_SKIP_LINE_PREFIXES = ("track", "browser", "#")
    VALID_INDEXINGS = ("0-based", "1-based")

    def __init__(self) -> None:
        """Constructor for the BedParser class."""

    def read_bed(
        self,
        filepath: str,
        names: list[str] | None = None,
        **kwargs: Any,
    ) -> pd.DataFrame | Any:
        """
        Load genomic intervals from a BED file.

        Ignores blank lines and UCSC header or comment lines that start with
        ``track``, ``browser``, or ``#``. Interval starts in the result are 0-based.

        For very large files, pass ``chunksize`` to read the intervals in batches
        instead of loading the whole table at once.

        :param filepath: Path to the BED file to load.
        :type filepath: str
        :param names: Column names for the BED fields. If ``None``, columns are numbered.
        :type names: list[str] | None
        :param kwargs: Extra options for ``pandas.read_csv`` (for example ``chunksize``).
        :return: Table of intervals, or a chunk iterator when ``chunksize`` is set.
        :rtype: pd.DataFrame | Any
        """
        with open(filepath) as handle:
            skiprows = [i for i, line in enumerate(handle) if not self._is_bed_data_line(line)]

        return pd.read_csv(
            filepath,
            sep="\t",
            header=None,
            names=names,
            skiprows=skiprows,
            **kwargs,
        )

    def _is_bed_data_line(self, line: str) -> bool:
        """
        Check whether a line is a BED interval row.

        Returns ``False`` for blank lines and for UCSC ``track``, ``browser``, or
        ``#`` comment lines.

        :param line: One line from a BED file.
        :type line: str
        :return: ``True`` if the line should be kept as interval data.
        :rtype: bool
        """
        return bool(line.strip()) and not line.startswith(self.BED_SKIP_LINE_PREFIXES)

    def write_bed(
        self,
        bed: pd.DataFrame,
        filepath: str,
        columns: list[str] | None = None,
    ) -> None:
        """
        Write genomic intervals to a BED file.

        The file has no header row. Starts must already be 0-based. If your
        coordinates come from FASTA headers or GFF/GTF (1-based), convert them
        with :meth:`convert_start` first.

        :param bed: Table of intervals to write.
        :type bed: pd.DataFrame
        :param filepath: Path of the BED file to create.
        :type filepath: str
        :param columns: Columns to write, in order. If ``None``, all columns are written.
        :type columns: list[str] | None
        :return: None
        """
        bed_to_write = bed[columns] if columns is not None else bed
        bed_to_write.to_csv(filepath, sep="\t", header=False, index=False)

    def convert_start(
        self,
        start_coordinates: Any,
        from_indexing: str,
        to_indexing: str,
    ) -> Any:
        """
        Convert interval start positions between 0-based and 1-based indexing.

        Use this when moving starts between BED (0-based) and formats such as
        GFF/GTF or region-generator FASTA headers (1-based). Ends are unchanged
        in the usual half-open BED and inclusive GFF conventions.

        :param start_coordinates: One start position or a series of starts.
        :type start_coordinates: Any
        :param from_indexing: Current indexing, ``0-based`` or ``1-based``.
        :type from_indexing: str
        :param to_indexing: Desired indexing, ``0-based`` or ``1-based``.
        :type to_indexing: str
        :return: Converted start(s), same type as ``start_coordinates``.
        :rtype: Any
        :raises ConfigurationError: If ``from_indexing`` or ``to_indexing`` is not supported.
        """
        if from_indexing not in self.VALID_INDEXINGS:
            raise ConfigurationError(
                f"Invalid from_indexing '{from_indexing}'. Expected one of {self.VALID_INDEXINGS}."
            )
        if to_indexing not in self.VALID_INDEXINGS:
            raise ConfigurationError(
                f"Invalid to_indexing '{to_indexing}'. Expected one of {self.VALID_INDEXINGS}."
            )
        if from_indexing == to_indexing:
            return start_coordinates
        if from_indexing == "0-based" and to_indexing == "1-based":
            return start_coordinates + 1
        return start_coordinates - 1
