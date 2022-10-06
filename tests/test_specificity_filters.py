import os

import pandas as pd

from oligo_designer_toolsuite.IO._data_parser import read_oligos_DB_tsv
from oligo_designer_toolsuite.oligo_specificity_filter._filter_blastn import (
    ProbeFilterBlastn,
)
from oligo_designer_toolsuite.oligo_specificity_filter._filter_bowtie import (
    ProbeFilterBowtie,
)

cwd = os.getcwd()

# Specify parameters
n_jobs = 1
ligation_region = 0
dir_output = cwd + "/output"
dir_annotations = cwd + "/data"
min_probes_per_gene = 2

# Reference transcriptome files for tests
file_transcriptome_fasta = dir_annotations + "/reference_sample.fna"
file_transcriptome_fasta2 = dir_annotations + "/reference_sample2.fna"

# Files containing probe info for tests
file_probe_info_match = dir_annotations + "/oligo_DB_match.tsv"
file_probe_info_no_match = dir_annotations + "/oligo_DB_no_match.tsv"

# blastn parameters
word_size = 10
percent_identity = 80
probe_length_min = 30
probe_length_max = 40
coverage = 50

# bowtie parameters
min_mismatches = 4
mismatch_region = 5

# Get probe info dictionary to use as input to filters
probe_info_dict_match = read_oligos_DB_tsv(file_probe_info_match)
probe_info_dict_no_match = read_oligos_DB_tsv(file_probe_info_no_match)

# Parameters to test ligation region argument
file_probe_info_ligation_match = dir_annotations + "/oligo_DB_ligation_match.tsv"
file_probe_info_ligation_nomatch = dir_annotations + "/oligo_DB_ligation_nomatch.tsv"

probe_info_dict_ligation_match = read_oligos_DB_tsv(file_probe_info_ligation_match)
probe_info_dict_ligation_no_match = read_oligos_DB_tsv(file_probe_info_ligation_nomatch)


def test_filter_bowtie_format():
    # Check that format of probe info dataframe is preserved when applying filter

    bowtie_filter = ProbeFilterBowtie(
        n_jobs,
        dir_output,
        dir_annotations,
        file_transcriptome_fasta,
        min_mismatches,
        mismatch_region,
        ligation_region,
    )

    bowtie_filter.apply(probe_info_dict_no_match)

    df_correct_format = pd.read_csv(file_probe_info_no_match)

    blast_sample_output = pd.read_csv(dir_output + "/probes_bowtie/probes_AGRN.txt")
    assert list(df_correct_format.columns) == list(blast_sample_output.columns)


def test_filter_bowtie_all_matches():
    # Check that bowtie filter filters out a sequence that is identified as a match for user-defined threshholds

    # Run blast filter
    bowtie_filter = ProbeFilterBowtie(
        n_jobs,
        dir_output,
        dir_annotations,
        file_transcriptome_fasta,
        min_mismatches,
        mismatch_region,
        ligation_region,
    )

    bowtie_filter.apply(probe_info_dict_match)

    # Check that gene of matching probe is added to file genes_with_insufficient_probes.txt
    blast_sample_output = pd.read_csv(
        dir_output + "/probes_bowtie/genes_with_insufficient_probes.txt",
        sep="\t",
        header=None,
    )
    assert not blast_sample_output.empty


def test_filter_blast_format():
    # Check that format of probe info dataframe is preserved when applying filter

    blast_filter = ProbeFilterBlastn(
        n_jobs,
        file_transcriptome_fasta,
        dir_output,
        dir_annotations,
        word_size,
        percent_identity,
        coverage,
        probe_length_min,
        probe_length_max,
        ligation_region,
    )

    blast_filter.apply(probe_info_dict_no_match)

    df_correct_format = pd.read_csv(file_probe_info_no_match)

    blast_sample_output = pd.read_csv(dir_output + "/probes_blast/probes_AGRN.txt")
    assert list(df_correct_format.columns) == list(blast_sample_output.columns)


def test_filter_blast_all_matches():
    # Check that blast filter filters out a sequence that is identified as a match for user-defined threshholds

    # Run blast filter
    blast_filter = ProbeFilterBlastn(
        n_jobs,
        file_transcriptome_fasta,
        dir_output,
        dir_annotations,
        word_size,
        percent_identity,
        coverage,
        probe_length_min,
        probe_length_max,
        ligation_region,
    )

    blast_filter.apply(probe_info_dict_match)

    # Check that gene of matching probe is added to file genes_with_insufficient_probes.txt
    blast_sample_output = pd.read_csv(
        dir_output + "/probes_blast/genes_with_insufficient_probes.txt",
        sep="\t",
        header=None,
    )
    assert not blast_sample_output.empty


def test_filter_ligation_bowtie_match():
    # Test that bowtie filter filters out probe where no mismatches are found in the ligation region

    bowtie_filter = ProbeFilterBowtie(
        n_jobs,
        dir_output,
        dir_annotations,
        file_transcriptome_fasta2,
        min_mismatches,
        mismatch_region,
        ligation_region=10,
    )

    bowtie_filter.apply(probe_info_dict_ligation_match)

    # Check that gene of matching probe is added to file genes_with_insufficient_probes.txt
    bowtie_sample_output = pd.read_csv(
        dir_output + "/probes_bowtie/genes_with_insufficient_probes.txt",
        sep="\t",
        header=None,
    )
    assert not bowtie_sample_output.empty


def test_filter_ligation_bowtie_no_match():
    # Test that bowtie filter keeps probe where atleast one mismatch is found in the ligation region

    bowtie_filter = ProbeFilterBowtie(
        n_jobs,
        dir_output,
        dir_annotations,
        file_transcriptome_fasta2,
        min_mismatches,
        mismatch_region,
        ligation_region=10,
    )

    bowtie_filter.apply(probe_info_dict_ligation_no_match)

    # Check that file genes_with_insufficient_probes.txt is empty
    bowtie_sample_output = pd.read_csv(
        dir_output + "/probes_bowtie/probes_WASH7P.txt",
        sep="\t",
        header=None,
    )
    assert not bowtie_sample_output.empty
