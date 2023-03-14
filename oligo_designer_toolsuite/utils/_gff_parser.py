############################################
# imports
############################################

import re
import gzip
import pickle
import warnings
import sys
import math

import pandas as pd
from collections import defaultdict

############################################
# GFF Parser Class
############################################


# class GffParser:
#     """
#     Parse a GFF3 or GTF file into a pandas DataFrame, attributes to columns.
#     Adapted from https://gist.github.com/slowkow/8101481
#     """

#     def __init__(self) -> None:
#         """Constructor method"""
#         self.GFF_HEADER = [
#             "seqid",
#             "source",
#             "type",
#             "start",
#             "end",
#             "score",
#             "strand",
#             "phase",
#         ]
#         self.R_SEMICOLON = re.compile(r"\s*;\s*")
#         self.R_COMMA = re.compile(r"\s*,\s*")
#         self.R_KEYVALUE = re.compile(r"(\s+|\s*=\s*)")

#         self.dataframe_gff = None

#     def read_gff(self, file: str, target_lines: int = None):
#         """Open an optionally gzipped GFF3/GTF file and return a pandas.DataFrame.

#         :param file: Filename of GFF3/GTF file.
#         :type file: str
#         :param target_lines: Read the first n lines or leave 'None' to read all lines, default: None
#         :type target_lines: int
#         :return: DataFrame with GFF§/GTF file content.
#         :rtype: pandas.DataFrame
#         """
#         # Each column is a list stored as a value in this dict.
#         result = defaultdict(list)

#         for i, line in enumerate(self._read_lines(file, target_lines)):
#             for key in line.keys():
#                 # This key has not been seen yet, so set it to None for all
#                 # previous lines.
#                 if key not in result:
#                     result[key] = [None] * i

#             # Ensure this row has some value for each column.
#             for key in result.keys():
#                 result[key].append(line.get(key, None))


#         result = pd.DataFrame(result)
#         cols = [
#             col
#             for col in result.columns
#             if col in self.GFF_HEADER or not result[col].isnull().all()
#         ]
#         dataframe = result[cols]

#         if self.dataframe_gff is not None:
#             warnings.warn(f"Overwriting existing gff dataframe!")
#         self.dataframe_gff = dataframe

#         return dataframe

#     def _read_lines(self, file: str, target_lines: int):
#         """Open an optionally gzipped GTF file and generate a dict for each line.

#         :param file: Filename of GFF3/GTF file.
#         :type file: str
#         :param target_lines: Read the first n lines or leave 'None' to read all lines, default: None
#         :type target_lines: int
#         :yield: Key - value pair of GFF3/GTF fileds or attributes
#         :rtype: dict
#         """
#         fn_open = gzip.open if file.endswith(".gz") else open
#         with fn_open(file) as fh:
#             for line_number, line in enumerate(fh):
#                 if target_lines is None or line_number <= target_lines:
#                     if line.startswith("#"):
#                         continue
#                     else:
#                         yield self._parse_fields(line)
#                 else:
#                     break

#     def _parse_fields(self, line: str):
#         """Parse a single GFF3/GTF line and return a dict.

#         :param line: Line of a GFF3/GTF file.
#         :type line: str
#         :return: Key - value pair of GFF3/GTF fileds or attributes.
#         :rtype: dict
#         """
#         result = {}

#         fields = line.rstrip().split("\t")

#         for i, col in enumerate(self.GFF_HEADER):
#             result[col] = self._get_value(fields[i])

#         # INFO field consists of "key1=value;key2=value;...".
#         infos = [x for x in re.split(self.R_SEMICOLON, fields[8]) if x.strip()]

#         for i, info in enumerate(infos, 1):
#             # It should be key="value".
#             try:
#                 key, _, value = re.split(self.R_KEYVALUE, info, 1)
#                 key = key.strip("\"'")
#             # But sometimes it is just "value".
#             except ValueError:
#                 key = "INFO{}".format(i)
#                 value = info
#             # Ignore the field if there is no value.
#             if value:
#                 result[key] = self._get_value(value)

#         return result

#     def _get_value(self, value: str):
#         """Return value for a field or attribute.

#         :param value: Field of a GFF3/GTF file.
#         :type value: str
#         :return: Value extracted from input field.
#         :rtype: str
#         """
#         if not value:
#             return None

#         # Strip double and single quotes.
#         value = value.strip("\"'")

#         # Return a list if the value has a comma.
#         if "," in value:
#             value = re.split(self.R_COMMA, value)
#         # These values are equivalent to None.
#         elif value in ["", ".", "NA"]:
#             return None

#         return value

#     def write_to_pickel(self, file_pickel):
#         """Store loaded GFF dataframe as pickel file.

#         :param file_pickel: File name of pickel file.
#         :type file_pickel: str
#         """
#         if self.dataframe_gff is None:
#             raise ValueError("No gff dataframe loaded!")

#         pickle.dump(self.dataframe_gff, file_pickel)

#     def read_from_pickel(self, file_pickel):
#         """Load GFF dataframe from pickel file

#         :param file_pickel:  File name of pickel file.
#         :type file_pickel: str
#         """
#         if self.dataframe_gff:
#             warnings.warn(f"Overwriting existing gff dataframe!")

#         self.dataframe_gff = pickle.load(open(file_pickel, "rb"))

#         return self.dataframe_gff


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
        
        self.dataframe_gff = None

    # Function to split GFF
    def _split_annotation(self, annotation_file, chunk_size, target_lines):
        csv_file = ''.join(annotation_file.split('.')[:-1]) + '.csv'
        extra_info_file = ''.join(annotation_file.split('.')[:-1]) + '.txt'
        
        finished = False
        lines_read = 0
        
        fn_open = gzip.open if annotation_file.endswith(".gz") else open
        with fn_open(annotation_file, 'r' ) as input_file:
            with open(csv_file, 'w') as out_csv:
                with open(extra_info_file, 'w') as out_extra_info:
                    while not finished and lines_read < target_lines:
                        content_chunck = ''
                        for _ in range(chunk_size):
                            lines_read += 1
                            if lines_read > target_lines:
                                break
                            try:
                                line = next(input_file)
                                if not line.startswith('#'):
                                    content_chunck += line
                            except:
                                finished = True
                        csv_content_chunck = re.sub('(.*)\t"(.*)"', r'\1', content_chunck)
                        extra_info_content_chunck = re.sub('(.*)\t"(.*)"', r'\2', content_chunck)

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
            # But sometimes it is just "value".
            except ValueError:
                key = "INFO{}".format(i)
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

        # Strip double and single quotes.
        value = value.strip("\"'")

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
        with open(info_file, 'r') as info_f:
            data = info_f.readlines()
            n_lines = len(data)
            for i in range(0, n_lines, chunk_size):
                data_chunk = data[i:i+chunk_size]
                info_dfs.append(self._info_to_df_chunk(data_chunk))
        return pd.concat(info_dfs)
    
    def read_gff(self, annotation_file, chunk_size=10000, target_lines=math.inf):
        csv_file, extra_info_file = self._split_annotation(annotation_file, chunk_size=chunk_size, target_lines=target_lines)
        
        info_df = self._info_to_df(extra_info_file, chunk_size=chunk_size)
        csv_df = pd.read_csv(csv_file, sep='\t', names=self.GFF_HEADER, header=None)
        
        csv_df.reset_index(inplace=True, drop=True)
        info_df.reset_index(inplace=True, drop=True)

        dataframe = pd.concat([csv_df, info_df], axis=1)
        
        if self.dataframe_gff is not None:
            warnings.warn(f"Overwriting existing gff dataframe!")
        self.dataframe_gff = dataframe

        return dataframe
    
    def write_to_pickel(self, file_pickel):
        """Store loaded GFF dataframe as pickel file.

        :param file_pickel: File name of pickel file.
        :type file_pickel: str
        """
        if self.dataframe_gff is None:
            raise ValueError("No gff dataframe loaded!")

        pickle.dump(self.dataframe_gff, file_pickel)

    def read_from_pickel(self, file_pickel):
        """Load GFF dataframe from pickel file

        :param file_pickel:  File name of pickel file.
        :type file_pickel: str
        """
        if self.dataframe_gff:
            warnings.warn(f"Overwriting existing gff dataframe!")

        self.dataframe_gff = pickle.load(open(file_pickel, "rb"))

        return self.dataframe_gff