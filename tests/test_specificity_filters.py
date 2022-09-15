import os

import pandas as pd

from oligo_designer_toolsuite.IO._data_parser import read_oligos_DB_tsv

# import sys
# sys.path.append("../oligo_designer_toolsuite")
from oligo_designer_toolsuite.oligo_specificity_filter._filter_blastn import (
    ProbeFilterBlastn,
)
from oligo_designer_toolsuite.oligo_specificity_filter._filter_bowtie import (
    ProbeFilterBowtie,
)

cwd = os.getcwd()

n_jobs = 10
ligation_region = 0
dir_output = cwd + "/output/"
dir_data = cwd + "/data/"
file_probe_info = (
    dir_data + "oligo_DB_unknown_unknown_Custom_release_unknown_gene_transcript.tsv"
)
min_probes_per_gene = 2
file_transcriptome_fasta = (
    dir_data
    + "reference_DB_unknown_unknown_Custom_release_unknown_genome_False_gene_transcript_True"
)
word_size = 10
percent_identity = 80
probe_length_min = 30
probe_length_max = 40
coverage = 50
min_mismatches = 4
mismatch_region = 5

probe_info_dict = read_oligos_DB_tsv(file_probe_info)


def test_filter_bowtie_format():

    # Run Bowtie filter
    bowtie_filter = ProbeFilterBowtie(
        n_jobs,
        dir_output,
        probe_length_min,
        probe_length_max,
        file_transcriptome_fasta,
        min_mismatches,
        mismatch_region,
    )

    bowtie_filter.apply(probe_info_dict)

    df_correct_format = pd.read_csv(file_probe_info)

    bowtie_sample_output = pd.read_csv(cwd + "/probes_bowtie/probes_AGRN.txt")
    assert (
        df_correct_format.columns().tolist() == bowtie_sample_output.columns().tolist()
    )

    # Give single sequence to filters
    # 1. Seq that is expected to be filtered out (Similar match)


def test_filter_blast_format():

    # Run Bowtie filter
    blast_filter = ProbeFilterBlastn(
        n_jobs,
        file_transcriptome_fasta,
        dir_output,
        word_size,
        percent_identity,
        coverage,
        probe_length_min,
        probe_length_max,
        ligation_region,
    )

    blast_filter.apply(probe_info_dict)

    df_correct_format = pd.read_csv(file_probe_info)

    blast_sample_output = pd.read_csv(dir_output + "/probes_blast/probes_AGRN.txt")
    assert (
        df_correct_format.columns().tolist() == blast_sample_output.columns().tolist()
    )


# TODO: Check if ligation filter works for all filters
