############################################
# imports
############################################

import os

from oligo_designer_toolsuite.database import OligoDatabase
from oligo_designer_toolsuite.oligo_specificity_filter import (
    Blastn,
    Bowtie,
    Bowtie2,
    BowtieSeedRegion,
    ExactMatches,
    LigationRegionCreation,
)

############################################
# Global Parameters
############################################

# Specify parameters
n_jobs = 1
ligation_region = 0
dir_annotations = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "data/specificity_filter"
)
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

############################################
# Tests
############################################


def test_filter_exact_matches():
    # check that exact matches filters out a doubled sequence from the db
    oligo_database = OligoDatabase()
    exact_matches = ExactMatches()
    oligo_database.load_database(file_oligo_info_exact_matches)
    filtered_oligo_info_dict_match = exact_matches.apply(
        oligo_database.database, None, n_jobs
    )

    assert (
        "WASH7P_1" not in filtered_oligo_info_dict_match["WASH7P"].keys()
    ), "A matching oligo has not been filtered from exact matches!"
    assert (
        "AGRN_1" not in filtered_oligo_info_dict_match["AGRN"].keys()
    ), "A matching oligo has not been filtered from exact mathces!"


def test_filter_bowtie_match(tmp_path):
    # Check that bowtie filter filters out a sequence that is identified as a match for user-defined threshholds

    oligo_database = OligoDatabase()
    bowtie_filter = Bowtie(
        tmp_path,
        num_mismatches,
        mismatch_region,
    )
    oligo_database.load_database(file_oligo_info_match)
    filtered_oligo_info_dict_match = bowtie_filter.apply(
        oligo_database.database, file_transcriptome_fasta, n_jobs
    )

    # check tha the oligo has been removed form the dataset
    assert (
        "WASH7P_1" not in filtered_oligo_info_dict_match["WASH7P"].keys()
    ), "A matching oligo has not been filtered from Bowtie!"


def test_filter_bowtie_no_match(tmp_path):
    # Check that bowtie does not filter filters out a sequence which is not a match
    oligo_database = OligoDatabase()
    bowtie_filter = Bowtie(
        tmp_path,
        num_mismatches,
        mismatch_region,
    )
    oligo_database.load_database(file_oligo_info_no_match)
    filtered_oligo_info_dict_match = bowtie_filter.apply(
        oligo_database.database, file_transcriptome_fasta, n_jobs
    )

    # check tha the oligo has been removed form the dataset
    assert (
        "AGRN_1" in filtered_oligo_info_dict_match["AGRN"].keys()
    ), "A non matching oligo has been filtered from Bowtie!"


def test_filter_bowtie2_match(tmp_path):
    # Check that bowtie filter filters out a sequence that is identified as a match for user-defined threshholds
    oligo_database = OligoDatabase()
    bowtie2_filter = Bowtie2(
        tmp_path,
    )
    oligo_database.load_database(file_oligo_info_match)
    filtered_oligo_info_dict_match = bowtie2_filter.apply(
        oligo_database.database, file_transcriptome_fasta, n_jobs
    )

    # check tha the oligo has been removed form the dataset
    assert (
        "WASH7P_1" not in filtered_oligo_info_dict_match["WASH7P"].keys()
    ), "A matching oligo has not been filtered from Bowtie!"


def test_filter_bowtie2_no_match(tmp_path):
    # Check that bowtie does not filter filters out a sequence which is not a match
    oligo_database = OligoDatabase()
    bowtie2_filter = Bowtie2(
        tmp_path,
    )
    oligo_database.load_database(file_oligo_info_no_match)
    filtered_oligo_info_dict_match = bowtie2_filter.apply(
        oligo_database.database, file_transcriptome_fasta, n_jobs
    )

    # check tha the oligo has been removed form the dataset
    assert (
        "AGRN_1" in filtered_oligo_info_dict_match["AGRN"].keys()
    ), "A non matching oligo has been filtered from Bowtie!"


def test_filter_blast_match(tmp_path):
    # Check that blast filter filters out a sequence that is identified as a match for user-defined threshholds

    # Run blast filter
    oligo_database = OligoDatabase()
    blast_filter = Blastn(
        tmp_path,
        coverage,
        word_size=word_size,
        perc_identity=percent_identity,
        strand=strand,
    )
    oligo_database.load_database(file_oligo_info_match)
    filtered_oligo_info_dict_match = blast_filter.apply(
        oligo_database.database, file_transcriptome_fasta, n_jobs
    )

    # check tha the oligo has been removed form the dataset
    assert (
        "WASH7P_1" not in filtered_oligo_info_dict_match["WASH7P"].keys()
    ), "A matching oligo has not been filtered from Blast!"


def test_filter_blast_no_match(tmp_path):
    # Check that blast does not filter filters out a sequence which is not a match
    oligo_database = OligoDatabase()
    blast_filter = Blastn(
        tmp_path,
        coverage,
        word_size=word_size,
        perc_identity=percent_identity,
        strand=strand,
    )
    oligo_database.load_database(file_oligo_info_no_match)
    filtered_oligo_info_dict_match = blast_filter.apply(
        oligo_database.database, file_transcriptome_fasta, n_jobs
    )

    # check tha the oligo has been removed form the dataset
    assert (
        "AGRN_1" in filtered_oligo_info_dict_match["AGRN"].keys()
    ), "A non matching oligo has been filtered from Bowtie!"


def test_seed_filter_match(tmp_path):
    oligo_database = OligoDatabase()
    ligation_seed_region = LigationRegionCreation(ligation_region_size=10)
    seed_region_filter = BowtieSeedRegion(tmp_path, ligation_seed_region)
    oligo_database.load_database(file_oligo_info_ligation_match)
    filtered_oligo_info_dict_ligation_match = seed_region_filter.apply(
        oligo_database.database, file_transcriptome_fasta_ligation, n_jobs
    )

    # check tha the oligo has been removed form the dataset
    assert (
        "WASH7P_1" not in filtered_oligo_info_dict_ligation_match["WASH7P"].keys()
    ), "A  matching oligo hasn't been filtered from Seed Region Filter!"


def test_seed_filter_no_match(tmp_path):
    oligo_database = OligoDatabase()
    ligation_seed_region = LigationRegionCreation(ligation_region_size=10)
    seed_region_filter = BowtieSeedRegion(tmp_path, ligation_seed_region)
    oligo_database.load_database(file_oligo_info_ligation_nomatch)
    filtered_oligo_info_dict_ligation_no_match = seed_region_filter.apply(
        oligo_database.database, file_transcriptome_fasta_ligation, n_jobs
    )
    # check tha the oligo has been removed form the dataset
    assert (
        "WASH7P_1" in filtered_oligo_info_dict_ligation_no_match["WASH7P"].keys()
    ), "A non matching oligo has been filtered from Seed Region Filter!"
