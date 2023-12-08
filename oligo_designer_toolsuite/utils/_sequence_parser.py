############################################
# imports
############################################

import os
import re
import gzip
import pickle
import warnings

import math
import pandas as pd

from Bio import SeqIO


############################################
# GFF Parser Class
############################################


class GffParser:
    def __init__(self) -> None:
        """Constructor method"""
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

    def check_gff_format(self, file):
        """File check. Is the file a GFF3/ or GTF file?

        :param file: Path to file.
        :type file: str
        :return: Returns True if file has GFF3 or GTF format.
        :rtype: bool
        """
        gtf = self.parse_annotation_from_gff(file, target_lines=100)
        return any(gtf)

    def parse_annotation_from_gff(
        self, annotation_file, file_pickel=None, chunk_size=10000, target_lines=math.inf
    ):
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

    def load_annotation_from_pickel(self, file_pickel):
        """Load GFF dataframe from pickel file

        :param file_pickel:  File name of pickel file.
        :type file_pickel: str
        """
        dataframe_gff = pickle.load(open(file_pickel, "rb"))

        return dataframe_gff

    # Function to split GFF
    def _split_annotation(self, annotation_file, chunk_size, target_lines):
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

    def _parse_fields(self, line: str):
        """Parse a single GFF3/GTF line and return a dict.

        :param line: Line of a GFF3/GTF file.
        :type line: str
        :return: Key - value pair of GFF3/GTF fileds or attributes.
        :rtype: dict
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
        """Return value for a field or attribute.

        :param value: Field of a GFF3/GTF file.
        :type value: str
        :return: Value extracted from input field.
        :rtype: str
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

    def _info_to_df_chunk(self, data_chunk):
        data_chunk = list(map(self._parse_fields, data_chunk))
        return pd.DataFrame(data_chunk)

    def _info_to_df(self, info_file, chunk_size):
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
        """Constructor method"""

    def check_fasta_format(self, file):
        """File check. Is the file a fasta file?

        :param file: Path to file.
        :type file: str
        :return: Returns True if file has fasta format.
        :rtype: bool
        """
        # taken from https://stackoverflow.com/questions/44293407/how-can-i-check-whether-a-given-file-is-fasta
        with open(file, "r") as handle:
            fasta = SeqIO.parse(handle, "fasta")
            return any(fasta)  # False when `fasta` is empty, i.e. wasn't a FASTA file

    def is_coordinate(self, entry):
        pattern = r"\S+:\S+-\S+\(.*\)"
        return bool(re.match(pattern, entry))

    def read_fasta_sequences(self, file_fasta_in, region_ids):
        # check if files exist and are in correct format
        if os.path.exists(file_fasta_in):
            if not self.check_fasta_format(file_fasta_in):
                raise ValueError("Fasta file has incorrect format!")
        else:
            raise ValueError("Fasta file does not exist!")
        # remove undesired regions already at the beginning
        if region_ids:
            fasta_sequences = []
            with open(file_fasta_in, "r") as handle:
                for entry in SeqIO.parse(handle, "fasta"):
                    region, _, _ = self.parse_fasta_header(entry.id, parse_additional_info=False)
                    if region in region_ids:
                        fasta_sequences.append(entry)
        else:
            with open(file_fasta_in, "r") as handle:
                fasta_sequences = list(SeqIO.parse(handle, "fasta"))

        return fasta_sequences

    def parse_fasta_header(self, header, parse_additional_info=True):
        """Helper function to parse the header of a sequence entry in a fasta file.

        :param header: Header of sequence entry, starting with '>'.
        :type header: string
        :return: Parsed header information, i.e. region, coordinates (chromosome, start, end, strand), and (optional) additional information
        :rtype: str, dict, str
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
                    coordinates.setdefault("chromosome", []).append(header_coordinate.split(":")[0])
                    coordinates.setdefault("start", []).append(
                        int(header_coordinate.split(":")[1].split("-")[0])
                    )
                    coordinates.setdefault("end", []).append(
                        int(header_coordinate.split(":")[1].split("-")[1].split("(")[0])
                    )
                    coordinates.setdefault("strand", []).append(header_coordinate.split("(")[1].split(")")[0])
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
