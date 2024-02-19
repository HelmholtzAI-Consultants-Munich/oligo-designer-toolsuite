############################################
# imports
############################################

import gzip
import math
import os
import pickle
import re

import pandas as pd
from Bio import SeqIO

############################################
# GFF Parser Class
############################################


class GffParser:
    """A class for parsing GFF (General Feature Format) files.

    GFF is a standard format for describing features of DNA, RNA, and protein sequences.
    This parser is designed to handle GFF files and provides methods for extracting information
    from the parsed data.
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

    def check_gff_format(self, file: str):
        """
        Checks the format of a GFF (General Feature Format) file.

        This method parses the GFF file up to the specified number of lines (target_lines) to determine if it
        contains any valid GFF entries.

        :param file: The path to the GFF file to be checked.
        :type file: str

        :raises ValueError: If the file is not a valid GFF file or does not exist.
        """

        def _check_gff_content(file):
            gtf = self.parse_annotation_from_gff(file, target_lines=100)
            return any(gtf)

        if os.path.exists(file):
            if not _check_gff_content(file):
                raise ValueError("GFF file has incorrect format!")
        else:
            raise ValueError("GFF file does not exist!")

    def parse_annotation_from_gff(
        self,
        annotation_file: str,
        file_pickel: str = None,
        chunk_size: int = 10000,
        target_lines: int = math.inf,
    ):
        """Parses annotation data from a GFF (General Feature Format) file.

        This method reads the GFF file, splits it into chunks for processing, and returns a DataFrame
        containing the parsed annotation data or the path to a pickle file with the parsed annotation data.

        :param annotation_file: The path to the GFF file to be parsed.
        :type annotation_file: str
        :param file_pickel: Optional. The path to save the pickled DataFrame. If None, the DataFrame is returned.
        :type file_pickel: Optional[str]
        :param chunk_size: The number of lines to process in each chunk.
        :type chunk_size: int
        :param target_lines: The maximum number of lines to read from the GFF file.
        :type target_lines: int
        :return: The parsed annotation data as a DataFrame or the path to the pickled DataFrame if 'file_pickel' is provided.
        :rtype: Union[pd.DataFrame, str]
        """
        csv_file, extra_info_file = self._split_annotation(
            annotation_file, chunk_size=chunk_size, target_lines=target_lines
        )

        info_df = self._info_to_df(extra_info_file, chunk_size=chunk_size)
        csv_df = pd.read_csv(csv_file, sep="\t", names=self.GFF_HEADER, header=None)

        csv_df.reset_index(inplace=True, drop=True)
        info_df.reset_index(inplace=True, drop=True)

        dataframe = pd.concat([csv_df, info_df], axis=1)

        os.remove(csv_file)
        os.remove(extra_info_file)

        if file_pickel:
            with open(file_pickel, "wb") as handle:
                pickle.dump(dataframe, handle)
            return file_pickel
        else:
            return dataframe

    def load_annotation_from_pickel(self, file_pickel: str):
        """
        Load annotation data from a pickled DataFrame.

        :param file_pickel: The path to the pickled DataFrame file.
        :type file_pickel: str
        :return: The loaded annotation data as a DataFrame.
        :rtype: pd.DataFrame
        """
        dataframe_gff = pickle.load(open(file_pickel, "rb"))

        return dataframe_gff

    def _split_annotation(
        self, annotation_file: str, chunk_size: int, target_lines: int
    ):
        """Split the GFF/GTF annotation file into a CSV file and an extra info file.

        :param annotation_file: The path to the GFF/GTF annotation file.
        :type annotation_file: str
        :param chunk_size: The number of lines to process in each chunk.
        :type chunk_size: int
        :param target_lines: The maximum number of lines to process.
        :type target_lines: int
        :return: The paths to the CSV file and the extra info file.
        :rtype: Tuple[str, str]
        """
        csv_file = ".".join(annotation_file.split(".")[:-1]) + ".csv"
        extra_info_file = ".".join(annotation_file.split(".")[:-1]) + ".txt"

        finished = False
        lines_read = 0

        fn_open = gzip.open if annotation_file.endswith(".gz") else open
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
                                    csv_content_chunck += (
                                        "\t".join(line.split("\t")[:8]) + "\n"
                                    )
                                    extra_info_content_chunck += "\t".join(
                                        line.split("\t")[8:]
                                    )
                            except:
                                finished = True

                        out_csv.write(csv_content_chunck)
                        out_extra_info.write(extra_info_content_chunck)

        return csv_file, extra_info_file

    def _parse_fields(self, line: str):
        """Parse the fields from a GFF/GTF annotation line.

        :param line: A line from the GFF/GTF annotation file.
        :type line: str
        :return: A dictionary containing parsed fields.
        :rtype: Dict[str, Any]
        """
        result = {}

        # INFO field consists of "key1=value;key2=value;...".
        infos = [x for x in re.split(self.R_SEMICOLON, line) if x.strip()]

        for i, info in enumerate(infos, 1):
            # It should be key="value".
            try:
                key, _, value = re.split(self.R_KEYVALUE, info, 1)
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

    def _get_value(self, value: str):
        """Get the parsed value from the GFF/GTF annotation field.

        :param value: The value of the GFF/GTF annotation field.
        :type value: str
        :return: Parsed value.
        :rtype: Union[None, str, List[str]]
        """
        if not value:
            return None

        # Strip double and single quotes and new lines.
        value = value.strip("\"'").strip("\n")

        # Return a list if the value has a comma.
        if "," in value:
            value = re.split(self.R_COMMA, value)
        # These values are equivalent to None.
        elif value in ["", ".", "NA"]:
            return None

        return value

    def _info_to_df_chunk(self, data_chunk: pd.DataFrame):
        """Convert a chunk of GFF/GTF annotation information to a DataFrame.

        :param data_chunk: Chunk of GFF/GTF annotation information.
        :type data_chunk: pd.DataFrame
        :return: DataFrame containing parsed information.
        :rtype: pd.DataFrame
        """
        data_chunk = list(map(self._parse_fields, data_chunk))
        return pd.DataFrame(data_chunk)

    def _info_to_df(self, info_file, chunk_size):
        """Convert GFF/GTF annotation information from a file to a DataFrame.

        :param info_file: Path to the file containing additional GFF/GTF annotation information.
        :type info_file: str
        :param chunk_size: Size of each chunk to read from the file.
        :type chunk_size: int
        :return: DataFrame containing parsed information.
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
    """A class for parsing FASTA files.

    This class provides methods for parsing FASTA files, extracting information
    from the headers, and handling sequence data.
    """

    def __init__(self) -> None:
        """Constructor for the FastaParser class."""

    def check_fasta_format(self, file: str):
        """Check the format of a FASTA file.

        This method verifies whether the given file is a valid FASTA file by attempting
        to parse its content using Biopython's SeqIO module.

        :param file: The path to the FASTA file to be checked.
        :type file: str

        :raises ValueError: If the file is not a valid FASTA file or does not exist.
        """

        def _check_fasta_content(file):
            """Check the content of a file to determine if it follows the FASTA format.

            This function reads the content of the provided file and uses Biopython's SeqIO
            module to attempt parsing it as a FASTA file. It returns True if the file contains
            at least one FASTA record, indicating a valid FASTA format. Otherwise, it returns False.

            :param file: The path to the file to be checked.
            :type file: str
            :return: True if the file is in FASTA format, False otherwise.
            :rtype: bool
            """
            # taken from https://stackoverflow.com/questions/44293407/how-can-i-check-whether-a-given-file-is-fasta
            with open(file, "r") as handle:
                fasta = SeqIO.parse(handle, "fasta")
                return any(
                    fasta
                )  # False when `fasta` is empty, i.e. wasn't a FASTA file

        if os.path.exists(file):
            if not _check_fasta_content(file):
                raise ValueError("Fasta file has incorrect format!")
            else:
                return True
        else:
            raise ValueError("Fasta file does not exist!")

    def is_coordinate(self, entry: str):
        """Check if the provided entry matches a specific coordinate pattern.

        This function uses a regular expression pattern to check if the entry follows
        the format of genomic coordinates, such as "chr:start-end(strand)". It returns
        True if the entry matches the pattern, indicating it represents genomic coordinates.

        :param entry: The entry to be checked for the coordinate pattern.
        :type entry: str
        :return: True if the entry is a genomic coordinate, False otherwise.
        :rtype: bool
        """
        pattern = r"\S+:\S+-\S+\(.*\)"
        return bool(re.match(pattern, entry))

    def get_fasta_regions(self, file_fasta_in: str):
        """Extract unique region identifiers from the headers of a FASTA file.

        This function parses the headers of a FASTA file to extract region identifiers.
        The region identifiers are then returned as a list, with duplicates removed.

        :param file_fasta_in: The input FASTA file.
        :type file_fasta_in: str
        :return: A list of unique region identifiers extracted from the FASTA file.
        :rtype: list[str]
        """
        region_ids = []
        with open(file_fasta_in, "r") as handle:
            for entry in SeqIO.parse(handle, "fasta"):
                region, _, _ = self.parse_fasta_header(
                    entry.id, parse_additional_info=False
                )
                region_ids.append(region)

        return list(set(region_ids))

    def read_fasta_sequences(self, file_fasta_in: str, region_ids: list[str] = None):
        """Read FASTA sequences from a file, optionally filtering by specified region identifiers.

        This function reads sequences from a FASTA file. If region_ids are provided, only the sequences
        corresponding to those regions will be included in the output. If no region_ids are provided,
        all sequences in the file will be returned.

        :param file_fasta_in: The input FASTA file.
        :type file_fasta_in: str
        :param region_ids: Optional list of region identifiers to filter the sequences.
        :type region_ids: list[str] or None
        :return: List of FASTA sequences from the file, filtered by region_ids if specified.
        :rtype: list[SeqRecord]
        """
        # remove undesired regions already at the beginning
        if region_ids:
            fasta_sequences = []
            with open(file_fasta_in, "r") as handle:
                for entry in SeqIO.parse(handle, "fasta"):
                    region, _, _ = self.parse_fasta_header(
                        entry.id, parse_additional_info=False
                    )
                    if region in region_ids:
                        fasta_sequences.append(entry)
        else:
            with open(file_fasta_in, "r") as handle:
                fasta_sequences = list(SeqIO.parse(handle, "fasta"))

        return fasta_sequences

    def parse_fasta_header(self, header: str, parse_additional_info: bool = True):
        """Parse information from a FASTA header.

        This function extracts region, additional information, and coordinates from a FASTA header.
        The header can contain region information, coordinate information, and additional
        information fields.

        :param header: The header string from a FASTA entry.
        :type header: str
        :param parse_additional_info: Flag to indicate whether to parse additional information.
                                    Default is True.
        :type parse_additional_info: bool
        :return: A tuple containing region, additional_info, and coordinates.
                - region: The region information extracted from the header.
                - additional_info: Additional information extracted from the header.
                - coordinates: Dictionary containing chromosome, start, end, and strand information.
        :rtype: tuple[str, dict, dict]
        """
        region = None
        additional_info = {}
        coordinates = {
            "chromosome": [None],
            "start": [None],
            "end": [None],
            "strand": [None],
        }

        for header_entry in header.split("::"):
            header_entry = header_entry.strip()
            if not region:
                region = header_entry
            elif self.is_coordinate(header_entry):
                header_coordinates = header_entry.split(";")
                coordinates = {}
                for header_coordinate in header_coordinates:
                    coordinates.setdefault("chromosome", []).append(
                        header_coordinate.split(":")[0]
                    )
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
                info_list = header_entry
                # the additional info field should be parsed, save information in dict
                if parse_additional_info:
                    if ";" in info_list:
                        info_list = info_list.split(";")

                        for infos in info_list:
                            key_values = infos.split(",")

                            for key_value in key_values:
                                key, value = key_value.split("=")

                                if key in additional_info:
                                    additional_info[key].append(value)
                                else:
                                    additional_info[key] = [value]

                    else:
                        key, value = info_list.split("=")
                        additional_info[key] = [value]
                else:
                    additional_info = info_list

        return region, additional_info, coordinates
