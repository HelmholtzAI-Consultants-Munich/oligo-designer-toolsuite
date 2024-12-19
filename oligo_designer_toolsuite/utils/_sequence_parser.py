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
        dataframe_gff = pickle.load(open(file_pickel, "rb"))

        return dataframe_gff

    def _split_annotation(self, annotation_file: str, chunk_size: int, target_lines: int) -> Tuple[str, str]:
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

    def _get_value(self, value: str) -> Union[None, str]:

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
        data_chunk = list(map(self._parse_fields, data_chunk))
        return pd.DataFrame(data_chunk)

    def _info_to_df(self, info_file: str, chunk_size: int) -> pd.DataFrame:
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

    def __init__(self) -> None:
        """Constructor for the FastaParser class."""

    def check_fasta_format(self, file: str) -> bool:
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
        pattern = r"\S+:\S+-\S+\(.*\)"
        return bool(re.match(pattern, entry))

    def get_fasta_regions(self, file_fasta_in: str) -> list:

        region_ids = []
        # use index instead of parse function for memory efficiency
        for idx in SeqIO.index(file_fasta_in, "fasta"):
            region, _, _ = self.parse_fasta_header(idx, parse_additional_info=False)
            region_ids.append(region)

        return list(set(region_ids))

    def read_fasta_sequences(self, file_fasta_in: str, region_ids: List[str] = None) -> list:
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
        files_in = check_if_list(files_in)
        file_out_mode = "w" if overwrite else "a"
        with open(file_out, file_out_mode) as out:
            for file_in in files_in:
                if os.path.isfile(file_in):
                    with open(file_in, "r") as in_f:
                        out.write(in_f.read())
