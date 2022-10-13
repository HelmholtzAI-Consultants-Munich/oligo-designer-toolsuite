import copy
import csv
import logging
import os
import shutil
from collections import OrderedDict
from os.path import exists

import numpy as np
import pandas as pd
import pybedtools
from BCBio import GFF
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from six import string_types
from six.moves import intern

# Adapted from https://github.com/openvax/gtfparse/blob/master/gtfparse/read_gtf.py


def read_gtf(
    filepath_or_buffer,
    expand_attribute_column=True,
    infer_biotype_column=False,
    column_converters={},
    usecols=None,
    features=None,
    chunksize=1024 * 1024,
):
    """
    Parse a GTF into a dictionary, mapping column names to sequences of values.

    :param filepath_or_buffer: Path to GTF file (may be gzip compressed) or buffer object such as StringIO
    :type filepath_or_buffer: str or buffer object
    :param expand_attribute_column: Replace strings of semi-colon separated key-value values in the
        'attribute' column with one column per distinct key, with a list of
        values for each row (using None for rows where key didn't occur).
    :type expand_attribute_column: bool
    :param infer_biotype_column:
        Due to the annoying ambiguity of the second GTF column across multiple
        Ensembl releases, figure out if an older GTF's source column is actually
        the gene_biotype or transcript_biotype.
    :type infer_biotype_column: bool
    :param column_converters:
        Dictionary mapping column names to conversion functions. Will replace
        empty strings with None and otherwise passes them to given conversion
        function.
    :type column_converters: dict, optional
    :param usecols:
        Restrict which columns are loaded to the give set. If None, then
        load all columns.
    :type usecols: list of str or None
    :param features:
        Drop rows which aren't one of the features in the supplied set
    :type features: set of str or None
    :param chunksize: Return the data in chunks of size chunksize
    :type chunksize: int
    """

    if isinstance(filepath_or_buffer, string_types) and not exists(filepath_or_buffer):
        raise ValueError("GTF file does not exist: %s" % filepath_or_buffer)

    if expand_attribute_column:
        result_df = parse_gtf_and_expand_attributes(
            filepath_or_buffer, chunksize=chunksize, restrict_attribute_columns=usecols
        )
    else:
        result_df = parse_gtf(result_df, features=features)

    for column_name, column_type in list(column_converters.items()):
        result_df[column_name] = [
            column_type(string_value) if len(string_value) > 0 else None
            for string_value in result_df[column_name]
        ]

    # Infer whether the values in the 'source' column of this GTF
    # are actually representing a biotype by checking for the most common
    # gene_biotype and transcript_biotype value 'protein_coding'
    if infer_biotype_column:
        unique_source_values = set(result_df["source"])
        if "protein_coding" in unique_source_values:
            column_names = set(result_df.columns)
            # Disambiguate between the two biotypes by checking if
            # gene_biotype is already present in another column. If it is,
            # the 2nd column is the transcript_biotype (otherwise, it's the
            # gene_biotype)
            if "gene_biotype" not in column_names:
                logging.info("Using column 'source' to replace missing 'gene_biotype'")
                result_df["gene_biotype"] = result_df["source"]
            if "transcript_biotype" not in column_names:
                logging.info(
                    "Using column 'source' to replace missing 'transcript_biotype'"
                )
                result_df["transcript_biotype"] = result_df["source"]

    if usecols is not None:
        column_names = set(result_df.columns)
        valid_columns = [c for c in usecols if c in column_names]
        result_df = result_df[valid_columns]

    return result_df


def parse_gtf_and_expand_attributes(
    filepath_or_buffer,
    chunksize=1024 * 1024,
    restrict_attribute_columns=None,
    features=None,
):
    """
    Parse lines into column->values dictionary and then expand
    the 'attribute' column into multiple columns. This expansion happens
    by replacing strings of semi-colon separated key-value values in the
    'attribute' column with one column per distinct key, with a list of
    values for each row (using None for rows where key didn't occur).

    :param filepath_or_buffer: Path to GTF file (may be gzip compressed) or buffer object such as StringIO
    :type filepath_or_buffer: str or buffer object
    :param chunksize: Return the data in chunks of size chunksize
    :type chunksize: int
    :param restrict_attribute_columns: If given, then only usese attribute columns.
    :type restrict_attribute_columns: list/set of str or None
    :param features: Ignore entries which don't correspond to one of the supplied features.
    :type features: set or None

    """
    result = parse_gtf(filepath_or_buffer, chunksize=chunksize, features=features)
    attribute_values = result["attribute"]
    del result["attribute"]
    for column_name, values in expand_attribute_strings(
        attribute_values, usecols=restrict_attribute_columns
    ).items():
        result[column_name] = values
    return result


def parse_gtf(
    filepath_or_buffer,
    chunksize=1024 * 1024,
    features=None,
    intern_columns=["seqname", "source", "strand", "frame"],
    fix_quotes_columns=["attribute"],
):
    """
    Parse gtf file

    :param filepath_or_buffer: Path to GTF file (may be gzip compressed) or buffer object such as StringIO
    :type filepath_or_buffer: str or buffer object
    :param chunksize: Return the data in chunks of size chunksize
    :type chunksize: int
    :param features: Ignore entries which don't correspond to one of the supplied features.
    :type features: set or None
    :param intern_columns: These columns are short strings which should be interned
    :type intern_columns:  list
    :param fix_quotes_columns: Most commonly the 'attribute' column which had broken quotes on some Ensembl release GTF files.
    :type fix_quotes_columns: list

    """

    if features is not None:
        features = set(features)

    dataframes = []

    def parse_frame(s):
        if s == ".":
            return 0
        else:
            return int(s)

    # GTF columns:
    # 1) seqname: str ("1", "X", "chrX", etc...)
    # 2) source : str
    #      Different versions of GTF use second column as of:
    #      (a) gene biotype
    #      (b) transcript biotype
    #      (c) the annotation source
    #      See: https://www.biostars.org/p/120306/#120321
    # 3) feature : str ("gene", "transcript", &c)
    # 4) start : int
    # 5) end : int
    # 6) score : float or "."
    # 7) strand : "+", "-", or "."
    # 8) frame : 0, 1, 2 or "."
    # 9) attribute : key-value pairs separated by semicolons
    # (see more complete description in docstring at top of file)

    REQUIRED_COLUMNS = [
        "seqname",
        "source",
        "feature",
        "start",
        "end",
        "score",
        "strand",
        "frame",
        "attribute",
    ]

    chunk_iterator = pd.read_csv(
        filepath_or_buffer,
        sep="\t",
        comment="#",
        names=REQUIRED_COLUMNS,
        skipinitialspace=True,
        skip_blank_lines=True,
        on_bad_lines="error",
        chunksize=chunksize,
        engine="c",
        dtype={
            "start": np.int64,
            "end": np.int64,
            "score": np.float32,
            "seqname": str,
        },
        na_values=".",
        converters={"frame": parse_frame},
    )
    dataframes = []
    try:
        for df in chunk_iterator:
            for intern_column in intern_columns:
                df[intern_column] = [intern(str(s)) for s in df[intern_column]]

            # compare feature strings after interning
            if features is not None:
                df = df[df["feature"].isin(features)]

            for fix_quotes_column in fix_quotes_columns:
                # Catch mistaken semicolons by replacing "xyz;" with "xyz"
                # Required to do this since the Ensembl GTF for Ensembl
                # release 78 has mistakes such as:
                #   gene_name = "PRAMEF6;" transcript_name = "PRAMEF6;-201"
                df[fix_quotes_column] = [
                    s.replace(';"', '"').replace(";-", "-")
                    for s in df[fix_quotes_column]
                ]
            dataframes.append(df)
    except Exception as e:
        raise ParsingError(str(e))
    df = pd.concat(dataframes)
    return df


class ParsingError(Exception):
    pass


def expand_attribute_strings(
    attribute_strings, quote_char='"', missing_value="", usecols=None
):
    """
    The last column of GTF has a variable number of key value pairs
    of the format: "key1 value1; key2 value2;"
    Parse these into a dictionary mapping each key onto a list of values,
    where the value is None for any row where the key was missing.

    :param attribute_strings: Values of the gtf attribute column.
    :type attribute_strings: list of str
    :param quote_char: Quote character to remove from values
    :type quote_char: str
    :param missing_value: If an attribute is missing from a row, give it this value.
    :type missing_value: any
    :param usecols: If not None, then only expand columns included in this set, otherwise use all columns.
    :type usecols: list of str or None

    :return: OrderedDict of column->value list mappings, in the order they appeared in the attribute strings.
    :rtype: OrderedDict

    """
    n = len(attribute_strings)

    extra_columns = {}
    column_order = []

    # While parsing millions of repeated strings (e.g. "gene_id" and "TP53"),
    # we can save a lot of memory by making sure there's only one string
    # object per unique string. The canonical way to do this is using
    # the 'intern' function. One problem is that Py2 won't let you intern
    # unicode objects, so to get around this we call intern(str(...)).
    #
    # It also turns out to be faster to check interned strings ourselves
    # using a local dictionary, hence the two dictionaries below
    # and pair of try/except blocks in the loop.
    column_interned_strings = {}
    value_interned_strings = {}

    for (i, attribute_string) in enumerate(attribute_strings):
        for kv in attribute_string.split(";"):
            # We're slicing the first two elements out of split() because
            # Ensembl release 79 added values like:
            #   transcript_support_level "1 (assigned to previous version 5)";
            # ...which gets mangled by splitting on spaces.
            parts = kv.strip().split(" ", 2)[:2]

            if len(parts) != 2:
                continue

            column_name, value = parts

            try:
                column_name = column_interned_strings[column_name]
            except KeyError:
                column_name = intern(str(column_name))
                column_interned_strings[column_name] = column_name

            if usecols is not None and column_name not in usecols:
                continue

            try:
                column = extra_columns[column_name]
            except KeyError:
                column = [missing_value] * n
                extra_columns[column_name] = column
                column_order.append(column_name)

            value = (
                value.replace(quote_char, "") if value.startswith(quote_char) else value
            )

            try:
                value = value_interned_strings[value]
            except KeyError:
                value = intern(str(value))
                value_interned_strings[value] = value

            # if an attribute is used repeatedly then
            # keep track of all its values in a list
            old_value = column[i]
            if old_value is missing_value:
                column[i] = value
            else:
                column[i] = "%s,%s" % (old_value, value)

    return OrderedDict(
        (column_name, extra_columns[column_name]) for column_name in column_order
    )


def get_sequence_from_annotation(
    file_bed, file_reference_fasta, file_fasta, split=False
):
    """
    Get sequence for regions annotated in file_bed using file_reference_fasta and save the output as a fasta file in file_fasta.

    :param file_bed: Path to bed file with annotated genomic regions.
    :type file_bed: string
    :param file_reference_fasta: Path to fasta file with reference sequence, e.g. transcriptome.
    :type file_reference_fasta: string
    :param file_fasta: Path to fasta file where retrieved sequences are written to.
    :type file_fasta: string
    :param split: Use -split option of bedtools getfasta, defaults to False
    :type split: bool
    """

    annotation = pybedtools.BedTool(file_bed)
    genome_sequence = pybedtools.BedTool(file_reference_fasta)

    annotation = annotation.sequence(fi=genome_sequence, s=True, name=True, split=split)
    annotation.save_seqs(file_fasta)


def check_gtf_format(filepath):
    with open(filepath) as in_handle:
        gtf = GFF.parse(in_handle, target_lines=1000)
        return any(gtf)


def check_fasta_format(file):
    # taken from https://stackoverflow.com/questions/44293407/how-can-i-check-whether-a-given-file-is-fasta
    with open(file, "r") as handle:
        fasta = SeqIO.parse(handle, "fasta")
        return any(fasta)  # False when `fasta` is empty, i.e. wasn't a FASTA file


def check_tsv_format(file):
    with open(file, "r") as tsv:
        read_tsv = csv.reader(tsv, delimiter="\t")
        return any(read_tsv)


def merge_fasta(files_fasta, file_merged_fasta):

    if files_fasta == []:
        raise ValueError("No fasta files provided for merge.")

    with open(file_merged_fasta, "wb") as handle_DB:
        for file in files_fasta:
            if os.path.exists(file):
                shutil.copyfileobj(open(file, "rb"), handle_DB)
                # shutil.rmtree(file)
                os.remove(file)
            else:
                raise ValueError(f"Fasta file {file} does not exist!")


def read_oligos_DB_gtf(file_oligos_DB_gtf, file_oligos_DB_fasta):
    """
    Create the oligo db dictionary from a gtf file.

    :param file_oligos_DB_gtf: Path to the file.
    :type file_oligos_DB_gtf: str
    :param file_oligos_DB_fasta: Path to the file.
    :type file_oligos_DB_fasta: str
    :return: oligos_DB
    :rtype: dict
    :raises ValueError: When the file given has the wrong format.
    :raises ValueError: When the path given does not exist.
    """

    if os.path.exists(file_oligos_DB_gtf):
        if not check_gtf_format(file_oligos_DB_gtf):
            raise ValueError("Database has incorrect format!")
    else:
        raise ValueError("Database file does not exist!")

    if os.path.exists(file_oligos_DB_fasta):
        if not check_fasta_format(file_oligos_DB_fasta):
            raise ValueError("Database has incorrect format!")
    else:
        raise ValueError("Database file does not exist!")

    oligos_df = read_gtf(file_oligos_DB_gtf)
    oligos_fasta = SeqIO.parse(file_oligos_DB_fasta, "fasta")
    oligos_DB = {}
    # compute the additional columns
    columns_fixed = [
        "probe_sequence",
        "transcript_id",
        "exon_id",
        "chromosome",
        "start",
        "end",
        "strand",
        "length",
    ]
    columns_df = list(oligos_df.columns)
    # remove the columns we don't need
    columns_df.remove("seqname")
    columns_df.remove("source")
    columns_df.remove("feature")
    columns_df.remove("score")
    columns_df.remove("frame")
    columns_df.remove("gene_id")
    additional_columns = [
        column for column in columns_df if (column not in columns_fixed)
    ]
    # crated the oligos db
    current_gene = ""
    current_probe = ""
    for index in range(len(oligos_df)):
        row = oligos_df.iloc[index]
        if row["gene_id"] != current_gene:
            current_gene = row["gene_id"]
            oligos_DB[current_gene] = {}
        if row["seqname"] != current_probe:
            current_probe = row["seqname"]
            oligos_DB[current_gene][current_probe] = {}
            oligo_fasta = next(oligos_fasta)
            probe_sequence = oligo_fasta.seq
            assert (
                oligo_fasta.id == current_probe
            )  # check that the fasta and gtf file are in sync
            oligos_DB[current_gene][current_probe] = {}
            # add all the values
            oligos_DB[current_gene][current_probe]["probe_sequence"] = probe_sequence
            oligos_DB[current_gene][current_probe]["transcript_id"] = [
                row["transcript_id"]
            ]
            oligos_DB[current_gene][current_probe]["exon_id"] = [row["exon_id"]]
            oligos_DB[current_gene][current_probe]["chromosome"] = row["chromosome"]
            oligos_DB[current_gene][current_probe]["start"] = [int(row["start"])]
            oligos_DB[current_gene][current_probe]["end"] = [int(row["end"])]
            oligos_DB[current_gene][current_probe]["strand"] = row["strand"]
            oligos_DB[current_gene][current_probe]["length"] = int(row["length"])
            for column in additional_columns:
                oligos_DB[current_gene][current_probe][column] = float(row[column])
        else:
            # append the values saved as a list
            oligos_DB[current_gene][current_probe]["transcript_id"].append(
                row["transcript_id"]
            )
            oligos_DB[current_gene][current_probe]["exon_id"].append(row["exon_id"])
            oligos_DB[current_gene][current_probe]["start"].append(row["start"])
            oligos_DB[current_gene][current_probe]["end"].append(row["end"])
    return oligos_DB


def read_oligos_DB_tsv(file_oligos_DB_tsv):
    """
    Reads a previously generated oligos DB and saves it in the ``oligos_DB`` attribute as a dictionary.
    The order of columns is :

    +----------+----------------+---------+---------------+---------+------------+-------+-----+--------+--------+------------------+
    +----------+----------------+---------+---------------+---------+------------+-------+-----+--------+--------+------------------+
    | probe_id | probe_sequence | gene_id | transcript_id | exon_id | chromosome | start | end | strand | length | additional feat. |
    +----------+----------------+---------+---------------+---------+------------+-------+-----+--------+--------+------------------+

    all the additional info computed by the filtering class.

    :param file_oligos_DB_tsv: path of the oligos_DB file
    :type file_oligos_DB_tsv: str
    :return: oligos_DB
    :rtype: dict
    """

    def parse_line_tsv(line, oligos_DB, current_gene, add_features):
        """Parses the lines of the tsv file of the oligos db and puts teh data in the dictionary.

        :param line: current line of the tsv file
        :type line: str
        :param current_gene: gene until wich we have created the dictionary
        :type current_gene: str
        :return: current gene we are after this iteration
        :rtype: str
        """
        line = line.split("\t")
        line[-1] = line[-1][0:-1]
        if line[2] != current_gene:
            current_gene = line[2]
            oligos_DB[current_gene] = {}
        probe_id = line[0]
        # what if we have duplicated sequences?
        oligos_DB[current_gene][probe_id] = {}
        oligos_DB[current_gene][probe_id]["probe_sequence"] = Seq(line[1])
        oligos_DB[current_gene][probe_id]["transcript_id"] = line[3].split(";")
        oligos_DB[current_gene][probe_id]["exon_id"] = line[4].split(";")
        oligos_DB[current_gene][probe_id]["chromosome"] = line[5]
        oligos_DB[current_gene][probe_id]["start"] = list(map(int, line[6].split(";")))
        oligos_DB[current_gene][probe_id]["end"] = list(map(int, line[7].split(";")))
        oligos_DB[current_gene][probe_id]["strand"] = line[8]
        oligos_DB[current_gene][probe_id]["length"] = int(line[9])
        # retrive the remaining features if they were computed
        if add_features:
            for i, column in enumerate(columns[10:]):
                oligos_DB[current_gene][probe_id][column] = float(line[i + 10])
        return current_gene, oligos_DB

    if os.path.exists(file_oligos_DB_tsv):
        if not check_tsv_format(file_oligos_DB_tsv):
            raise ValueError("Database has incorrect format!")
    else:
        raise ValueError("Database file does not exist!")

    oligos_DB = {}
    handle_probe = open(file_oligos_DB_tsv, "r")
    # read the header
    line = handle_probe.readline()
    columns = line.split("\t")
    columns[-1] = columns[-1][0:-1]  # delete \n in the last word
    add_features = len(columns) > 10
    # read the rest of the file
    current_gene = ""
    for line in handle_probe:
        current_gene, oligos_DB = parse_line_tsv(
            line, oligos_DB, current_gene, add_features
        )

    handle_probe.close()
    return oligos_DB


def write_oligos_DB_gtf(oligos_DB, file_oligos_DB_gtf, file_oligos_DB_fasta):
    """
    Writes the data structure oligos_DB in a gtf file in the ``file_oligos_DB_gtf`` path.
    The additional features are written in the 9th column and the sequence of the probes is written on a separate fasta fila
    with heading the probe_id.

    :param oligos_DB: oligos_DB dictionary
    :type oligos_DB: dict
    :param file_oligos_DB_gtf: path of the gtf file
    :type file_oligos_DB_gtf: str
    :param file_oligos_DB_fasta: path of the fasta file
    :type file_oligos_DB_fasta: str
    """
    with open(file_oligos_DB_gtf, "w") as handle_gtf:
        with open(file_oligos_DB_fasta, "w") as handle_fasta:
            source = "oligo-designer-toolsuite"
            feature = "Oligonucleotide"
            score = "."
            frame = "."
            output_fasta = []
            for gene_id, probe in oligos_DB.items():
                for probe_id, probe_attributes in probe.items():
                    output_fasta.append(
                        SeqRecord(probe_attributes["probe_sequence"], probe_id, "", "")
                    )  # write the sequence in the fasta file
                    for i in range(len(probe_attributes["start"])):
                        # write the annotation file
                        start = str(probe_attributes["start"][i])
                        end = str(probe_attributes["end"][i])
                        strand = probe_attributes["strand"]
                        output = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t".format(
                            probe_id,
                            source,
                            feature,
                            start,
                            end,
                            score,
                            strand,
                            frame,
                        )
                        # add all the oher features
                        output += f'gene_id "{gene_id}"; '
                        probe_attributes_copy = copy.deepcopy(probe_attributes)
                        # delete already written features
                        del probe_attributes_copy["start"]
                        del probe_attributes_copy["end"]
                        del probe_attributes_copy["strand"]
                        del probe_attributes_copy["probe_sequence"]
                        for key, value in probe_attributes_copy.items():
                            if type(value) == list:
                                output += f'{key} "{value[i]}"; '
                            else:
                                output += f'{key} "{value}"; '
                        output += "\n"
                        handle_gtf.write(output)
            SeqIO.write(output_fasta, handle_fasta, "fasta")


def write_oligos_DB_tsv(oligos_DB, file_oligos_DB_tsv):
    """
    Writes the data structure self.oligos_DB in a tsv file in the ``file_oligos_DB_tsv`` path.
    The order of columns is:

    +----------+----------------+---------+---------------+---------+------------+-------+-----+--------+--------+------------------+
    +----------+----------------+---------+---------------+---------+------------+-------+-----+--------+--------+------------------+
    | probe_id | probe_sequence | gene_id | transcript_id | exon_id | chromosome | start | end | strand | length | additional feat. |
    +----------+----------------+---------+---------------+---------+------------+-------+-----+--------+--------+------------------+

    Where all the additional info computed by the filtering class.

    :param oligos_DB: oligos_DB dictionary
    :type oligos_DB: dict
    :param file_oligos_DB_tsv: path of the tsv file
    :type file_oligos_DB_tsv: str
    """

    with open(file_oligos_DB_tsv, "w") as handle_probe:
        columns = [
            "probe_id",
            "probe_sequence",
            "gene_id",
            "transcript_id",
            "exon_id",
            "chromosome",
            "start",
            "end",
            "strand",
            "length",
        ]  # keep fixed the structure for these coulums
        # find all the other names of the columns (depend on the filters applied)
        genes = list(oligos_DB.keys())
        i = 0
        while oligos_DB[genes[i]] == {}:
            i += 1
        tmp = oligos_DB[genes[i]]
        tmp = list(list(tmp.values())[0].keys())
        additional_columns = [column for column in tmp if (column not in columns)]
        columns.extend(additional_columns)
        handle_probe.write("\t".join(columns) + "\n")

        for gene_id, probe in oligos_DB.items():
            for probe_id, probe_attributes in probe.items():
                # write the basic information information we compute for each probe
                output = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
                    probe_id,
                    probe_attributes["probe_sequence"],
                    gene_id,
                    ";".join(probe_attributes["transcript_id"]),
                    ";".join(probe_attributes["exon_id"]),
                    probe_attributes["chromosome"],
                    ";".join(str(s) for s in probe_attributes["start"]),
                    ";".join(str(e) for e in probe_attributes["end"]),
                    probe_attributes["strand"],
                    probe_attributes["length"],
                )
                # if we computed additional features we write also those
                if len(additional_columns) > 0:
                    for column in additional_columns:
                        output += "\t{}".format(probe_attributes[column])
                    # \n at the end f the string
                output += "\n"
                handle_probe.write(output)
