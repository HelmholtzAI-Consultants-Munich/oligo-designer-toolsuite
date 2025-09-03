############################################
# imports
############################################

import gzip
import math
import os
import pickle
import re
from typing import List, Union, Tuple

import pandas as pd
from Bio import SeqIO

from oligo_designer_toolsuite._constants import (
    SEPARATOR_FASTA_HEADER_FIELDS,
    SEPARATOR_FASTA_HEADER_FIELDS_LIST,
    SEPARATOR_FASTA_HEADER_FIELDS_LIST_ITEMS,
)

from ._checkers_and_helpers import check_if_list

############################################
# GFF Parser Class
############################################


class GffParser:
    """
    A parser class for handling GFF (General Feature Format) files.

    GFF is a standard format for describing features of DNA, RNA, and protein sequences.
    This class provides utilities to parse and process GFF files.
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
        Checks if the provided GFF file is in the correct format by attempting to parse its content.

        :param file: The path to the GFF file to be checked.
        :type file: str
        :raises ValueError: If the file does not exist or if the GFF format is incorrect.
        :return: True if the file exists and is in the correct format, otherwise raises an error.
        :rtype: bool
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
        file_pickle: str = None,
        chunk_size: int = 10000,
        target_lines: int = math.inf,
    ) -> Union[str, pd.DataFrame]:
        """
        Parses the GFF annotation file and converts it into a DataFrame. Optionally, saves the DataFrame to a pickle file.

        :param annotation_file: The path to the GFF annotation file to be parsed.
        :type annotation_file: str
        :param file_pickle: The path to save the resulting DataFrame as a pickle file (optional).
        :type file_pickle: str, optional
        :param chunk_size: The number of lines to process at a time from the GFF file.
        :type chunk_size: int
        :param target_lines: The number of lines to parse before stopping (useful for sampling).
        :type target_lines: int
        :return: The parsed annotation as a DataFrame or the path to the pickle file if specified.
        :rtype: Union[str, pd.DataFrame]
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

        if file_pickle:
            with open(file_pickle, "wb") as handle:
                pickle.dump(dataframe, handle)
            return file_pickle
        else:
            return dataframe

    def load_annotation_from_pickle(self, file_pickel: str) -> pd.DataFrame:
        """
        Loads a GFF annotation DataFrame from a pickle file.

        :param file_pickel: The path to the pickle file containing the GFF annotation DataFrame.
        :type file_pickel: str
        :return: The loaded DataFrame containing the GFF annotation.
        :rtype: pd.DataFrame
        """
        dataframe_gff = pickle.load(open(file_pickel, "rb"))

        return dataframe_gff

    def _split_annotation(self, annotation_file: str, chunk_size: int, target_lines: int) -> Tuple[str, str]:
        """
        Splits the GFF annotation file into two separate files: one containing the standard GFF columns and the other containing the extra information.

        :param annotation_file: The path to the GFF annotation file.
        :type annotation_file: str
        :param chunk_size: The number of lines to read at a time from the annotation file.
        :type chunk_size: int
        :param target_lines: The maximum number of lines to process from the annotation file.
        :type target_lines: int
        :return: A tuple containing the paths to the CSV file with the main GFF content and the text file with extra information.
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
                                    csv_content_chunck += "\t".join(line.split("\t")[:8]) + "\n"
                                    extra_info_content_chunck += "\t".join(line.split("\t")[8:])
                            except:
                                finished = True

                        out_csv.write(csv_content_chunck)
                        out_extra_info.write(extra_info_content_chunck)

        return csv_file, extra_info_file

    def _parse_fields(self, line: str) -> dict:
        """
        Parse the attributes from a GFF line.

        :param line: The line of text containing attribute information from a GFF file.
        :type line: str
        :return: A dictionary where keys are attribute names and values are the corresponding values from the GFF line.
        :rtype: dict
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

    def _get_value(self, value: str) -> Union[None, str]:
        """
        Process the value extracted from a GFF line, handling special cases like lists, empty values, and standardization.

        :param value: The value from an attribute in a GFF line.
        :type value: str
        :return: Processed value, which may be None, a string, or a list of strings.
        :rtype: Union[None, str]
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

    def _info_to_df_chunk(self, data_chunk: pd.DataFrame) -> pd.DataFrame:
        """
        Convert a chunk of data into a DataFrame by parsing its fields.

        :param data_chunk: A chunk of data read from the GFF info file.
        :type data_chunk: pd.DataFrame
        :return: A DataFrame with parsed fields.
        :rtype: pd.DataFrame
        """
        data_chunk = list(map(self._parse_fields, data_chunk))
        return pd.DataFrame(data_chunk)

    def _info_to_df(self, info_file: str, chunk_size: int) -> pd.DataFrame:
        """
        Read an info file in chunks and convert it into a DataFrame.

        :param info_file: The path to the info file containing GFF attribute data.
        :type info_file: str
        :param chunk_size: The number of lines to process in each chunk.
        :type chunk_size: int
        :return: A DataFrame combining all chunks.
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
    A parser class for handling FASTA files.

    The FastaParser class provides methods for processing and manipulating sequences within FASTA files.
    It is designed to support various sequence-related tasks, such as reading sequences, parsing sequence headers,
    and performing operations on the parsed data.
    """

    def __init__(self) -> None:
        """Constructor for the FastaParser class."""

    def check_fasta_format(self, file: str) -> bool:
        """
        Validates whether the given file is in correct FASTA format.

        This function checks if a file exists and verifies whether it is in proper FASTA format by inspecting its content.

        :param file: The path to the file to be checked.
        :type file: str
        :return: True if the file is a valid FASTA file, otherwise raises a ValueError.
        :rtype: bool
        """

        def _check_fasta_content(file) -> bool:
            fasta = SeqIO.index(file, "fasta")
            return any(fasta)  # False when `fasta` is empty, i.e. wasn't a FASTA file

        if os.path.exists(file):
            if not _check_fasta_content(file):
                raise ValueError("Fasta file has incorrect format!")
            else:
                return True
        else:
            raise ValueError("Fasta file does not exist!")

    def is_coordinate(self, entry: str) -> bool:
        """
        Checks if a given string is a valid coordinate format.

        :param entry: The string to be checked for coordinate format.
        :type entry: str
        :return: True if the string matches the coordinate pattern, otherwise False.
        :rtype: bool
        """
        pattern = r"\S+:\S+-\S+\(.*\)"
        return bool(re.match(pattern, entry))

    def get_fasta_regions(self, file_fasta_in: str) -> list:
        """
        Extracts unique region identifiers from a FASTA file.

        :param file_fasta_in: The path to the input FASTA file.
        :type file_fasta_in: str
        :return: A list of unique region identifiers found in the FASTA file.
        :rtype: list
        """
        region_ids = []
        # use index instead of parse function for memory efficiency
        for idx in SeqIO.index(file_fasta_in, "fasta"):
            region, _, _ = self.parse_fasta_header(idx, parse_additional_info=False)
            region_ids.append(region)

        return list(set(region_ids))

    def read_fasta_sequences(self, file_fasta_in: str, region_ids: List[str] = None) -> list:
        """
        Reads sequences from a FASTA file, optionally filtering by specific region identifiers.

        This function reads sequences from the specified FASTA file. If a list of region IDs is provided,
        only the sequences corresponding to those regions are extracted.

        :param file_fasta_in: The path to the input FASTA file.
        :type file_fasta_in: str
        :param region_ids: List of region IDs to process. If None, all regions in the OligoDatabase are processed, defaults to None.
        :type region_ids: Union[str, List[str]], optional
        :return: A list of sequences from the FASTA file, filtered by region if specified.
        :rtype: list
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
    ) -> Tuple[str, dict, dict]:
        """
        Parses the header of a FASTA sequence to extract region, coordinates, and additional information.

        This function processes the header of a FASTA sequence to extract the region name,
        optional genomic coordinates, and additional key-value pair information.

        :param header: The header string from a FASTA sequence.
        :type header: str
        :param parse_coordinates: Whether to parse genomic coordinates from the header, defaults to True.
        :type parse_coordinates: bool, optional
        :param parse_additional_info: Whether to parse additional information from the header, defaults to True.
        :type parse_additional_info: bool, optional
        :return: A tuple containing the region name, additional info dictionary, and coordinates dictionary.
        :rtype: Tuple[str, dict, dict]
        """
        region = None
        additional_info = {}
        coordinates = {
            "chromosome": [None],
            "start": [None],
            "end": [None],
            "strand": [None],
        }

        for header_entry in header.split(SEPARATOR_FASTA_HEADER_FIELDS):
            header_entry = header_entry.strip()
            if not region:
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
                info_list = header_entry
                # the additional info field should be parsed, save information in dict
                if parse_additional_info:
                    if SEPARATOR_FASTA_HEADER_FIELDS_LIST in info_list:
                        info_list = info_list.split(SEPARATOR_FASTA_HEADER_FIELDS_LIST)

                        for infos in info_list:
                            key_values = infos.split(SEPARATOR_FASTA_HEADER_FIELDS_LIST_ITEMS)

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

    def merge_fasta_files(self, files_in: list, file_out: str, overwrite: bool = False) -> None:
        """
        Merges multiple FASTA files into a single output file.

        This method takes a list of input FASTA files and merges their contents into a single output file.
        The output file can either overwrite any existing file or append to it based on the `overwrite` parameter.

        :param files_in: A list of input FASTA files to be merged.
        :type files_in: list
        :param file_out: The path to the output FASTA file.
        :type file_out: str
        :param overwrite: Whether to overwrite the output file if it exists, defaults to False.
        :type overwrite: bool
        """
        files_in = check_if_list(files_in)
        file_out_mode = "w" if overwrite else "a"
        with open(file_out, file_out_mode) as out:
            for file_in in files_in:
                if os.path.isfile(file_in):
                    with open(file_in, "r") as in_f:
                        out.write(in_f.read())
