import os
import shutil

import pytest

from oligo_designer_toolsuite.oligo_specificity_filter import (
    Blastn,
    Bowtie,
    Bowtie2,
    BowtieSeedRegion,
    ExactMatches,
    LigationRegionCreation,
)
from oligo_designer_toolsuite.utils import read_oligos_DB_tsv

cwd = os.getcwd()

# Specify parameters
n_jobs = 1
ligation_region = 0
dir_output = cwd + "/tests/output"
dir_annotations = cwd + "/tests/data"
min_oligos_per_gene = 2

# Reference transcriptome files for tests
file_transcriptome_fasta = dir_annotations + "/reference_sample.fna"
file_transcriptome_fasta_ligation = dir_annotations + "/reference_sample_ligation.fna"

# Files containing oligo info for tests
file_oligo_info_match = dir_annotations + "/oligo_DB_match.tsv"
file_oligo_info_no_match = dir_annotations + "/oligo_DB_no_match.tsv"
file_oligo_info_exact_matches = dir_annotations + "/oligo_DB_exact_matches.tsv"

# blastn parameters
word_size = 10
percent_identity = 80
oligo_length_min = 30
oligo_length_max = 40
coverage = 50
strand = "plus"

# bowtie parameters
num_mismatches = 3
mismatch_region = 5


# Parameters to test ligation region argument
file_oligo_info_ligation_match = dir_annotations + "/oligo_DB_ligation_match.tsv"
file_oligo_info_ligation_nomatch = dir_annotations + "/oligo_DB_ligation_no_match.tsv"


@pytest.fixture(autouse=True)
def run_around_tess():
    yield
    shutil.rmtree(dir_output)


def test_filter_exact_matches():
    # check that exact matches filters out a doubled sequence from the db
    exact_matches = ExactMatches(dir_output)
    oligo_info_dict_exact_matches = read_oligos_DB_tsv(file_oligo_info_exact_matches)
    filtered_oligo_info_dict_match = exact_matches.apply(
        oligo_info_dict_exact_matches, file_transcriptome_fasta, n_jobs
    )

    assert (
        "WASH7P_1" not in filtered_oligo_info_dict_match["WASH7P"].keys()
    ), "A matching oligo has not been filtered from exact matches!"
    assert (
        "AGRN_1" not in filtered_oligo_info_dict_match["AGRN"].keys()
    ), "A matching oligo has not been filtered from exact mathces!"


def test_filter_bowtie_match():
    # Check that bowtie filter filters out a sequence that is identified as a match for user-defined threshholds

    bowtie_filter = Bowtie(
        dir_output,
        num_mismatches,
        mismatch_region,
    )
    oligo_info_dict_match = read_oligos_DB_tsv(file_oligo_info_match)
    filtered_oligo_info_dict_match = bowtie_filter.apply(
        oligo_info_dict_match, file_transcriptome_fasta, n_jobs
    )

    # check tha the oligo has been removed form the dataset
    assert (
        "WASH7P_1" not in filtered_oligo_info_dict_match["WASH7P"].keys()
    ), "A matching oligo has not been filtered from Bowtie!"


def test_filter_bowtie_no_match():
    # Check that bowtie does not filter filters out a sequence which is not a match
    bowtie_filter = Bowtie(
        dir_output,
        num_mismatches,
        mismatch_region,
    )
    oligo_info_dict_no_match = read_oligos_DB_tsv(file_oligo_info_no_match)
    filtered_oligo_info_dict_match = bowtie_filter.apply(
        oligo_info_dict_no_match, file_transcriptome_fasta, n_jobs
    )

    # check tha the oligo has been removed form the dataset
    assert (
        "AGRN_1" in filtered_oligo_info_dict_match["AGRN"].keys()
    ), "A non matching oligo has been filtered from Bowtie!"


def test_filter_bowtie2_match():
    # Check that bowtie filter filters out a sequence that is identified as a match for user-defined threshholds

    bowtie2_filter = Bowtie2(
        dir_output,
    )
    oligo_info_dict_match = read_oligos_DB_tsv(file_oligo_info_match)
    filtered_oligo_info_dict_match = bowtie2_filter.apply(
        oligo_info_dict_match, file_transcriptome_fasta, n_jobs
    )

    # check tha the oligo has been removed form the dataset
    assert (
        "WASH7P_1" not in filtered_oligo_info_dict_match["WASH7P"].keys()
    ), "A matching oligo has not been filtered from Bowtie!"


def test_filter_bowtie2_no_match():
    # Check that bowtie does not filter filters out a sequence which is not a match
    bowtie2_filter = Bowtie2(
        dir_output,
    )
    oligo_info_dict_no_match = read_oligos_DB_tsv(file_oligo_info_no_match)
    filtered_oligo_info_dict_match = bowtie2_filter.apply(
        oligo_info_dict_no_match, file_transcriptome_fasta, n_jobs
    )

    # check tha the oligo has been removed form the dataset
    assert (
        "AGRN_1" in filtered_oligo_info_dict_match["AGRN"].keys()
    ), "A non matching oligo has been filtered from Bowtie!"


def test_filter_blast_match():
    # Check that blast filter filters out a sequence that is identified as a match for user-defined threshholds

    # Run blast filter
    blast_filter = Blastn(dir_output, word_size, percent_identity, coverage, strand)
    oligo_info_dict_match = read_oligos_DB_tsv(file_oligo_info_match)
    filtered_oligo_info_dict_match = blast_filter.apply(
        oligo_info_dict_match, file_transcriptome_fasta, n_jobs
    )

    # check tha the oligo has been removed form the dataset
    assert (
        "WASH7P_1" not in filtered_oligo_info_dict_match["WASH7P"].keys()
    ), "A matching oligo has not been filtered from Blast!"


def test_filter_blast_no_match():
    # Check that blast does not filter filters out a sequence which is not a match
    blast_filter = Blastn(dir_output, word_size, percent_identity, coverage, strand)
    oligo_info_dict_no_match = read_oligos_DB_tsv(file_oligo_info_no_match)
    filtered_oligo_info_dict_match = blast_filter.apply(
        oligo_info_dict_no_match, file_transcriptome_fasta, n_jobs
    )

    # check tha the oligo has been removed form the dataset
    assert (
        "AGRN_1" in filtered_oligo_info_dict_match["AGRN"].keys()
    ), "A non matching oligo has been filtered from Bowtie!"


def test_seed_filter_match():
    ligation_seed_region = LigationRegionCreation(ligation_region_size=10)
    seed_region_filter = BowtieSeedRegion(dir_output, ligation_seed_region)
    oligo_info_dict_ligation_match = read_oligos_DB_tsv(file_oligo_info_ligation_match)
    filtered_oligo_info_dict_ligation_match = seed_region_filter.apply(
        oligo_info_dict_ligation_match, file_transcriptome_fasta_ligation, n_jobs
    )

    # check tha the oligo has been removed form the dataset
    assert (
        "WASH7P_1" not in filtered_oligo_info_dict_ligation_match["WASH7P"].keys()
    ), "A  matching oligo hasn't been filtered from Seed Region Filter!"


def test_seed_filter_no_match():
    ligation_seed_region = LigationRegionCreation(ligation_region_size=10)
    seed_region_filter = BowtieSeedRegion(dir_output, ligation_seed_region)
    oligo_info_dict_ligation_no_match = read_oligos_DB_tsv(
        file_oligo_info_ligation_nomatch
    )
    filtered_oligo_info_dict_ligation_no_match = seed_region_filter.apply(
        oligo_info_dict_ligation_no_match, file_transcriptome_fasta_ligation, n_jobs
    )
    # check tha the oligo has been removed form the dataset
    assert (
        "WASH7P_1" in filtered_oligo_info_dict_ligation_no_match["WASH7P"].keys()
    ), "A non matching oligo has been filtered from Seed Region Filter!"
