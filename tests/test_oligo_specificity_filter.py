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
    CrossHybridizationFilter,
    ExactMatches,
    LigationRegionCreation,
)
from oligo_designer_toolsuite.oligo_specificity_filter.cross_hybridization_policies import (
    remove_by_bigger_region_policy,
    remove_by_degree_policy,
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
file_oligo_info_cross_hybridization = (
    dir_annotations + "/oligo_DB_cross_hybridization.tsv"
)

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

# Parameters Cross-hybridization
expected_oligos_bigger_region = {
    "region_1": {
        "region_1::oligo_7",
        "region_1::oligo_5",
        "region_1::oligo_6",
        "region_1::oligo_8",
        "region_1::oligo_4",
    },
    "region_2": {
        "region_2::oligo_3",
        "region_2::oligo_2",
        "region_2::oligo_6",
        "region_2::oligo_5",
        "region_2::oligo_4",
    },
    "region_3": {
        "region_3::oligo_1",
        "region_3::oligo_4",
        "region_3::oligo_3",
        "region_3::oligo_2",
        "region_3::oligo_5",
    },
}

expected_oligos_degree = {
    "region_1": {
        "region_1::oligo_1",
        "region_1::oligo_4",
        "region_1::oligo_5",
        "region_1::oligo_6",
        "region_1::oligo_7",
    },
    "region_2": {
        "region_2::oligo_1",
        "region_2::oligo_2",
        "region_2::oligo_3",
        "region_2::oligo_4",
        "region_2::oligo_5",
        "region_2::oligo_6",
        "region_2::oligo_7",
    },
    "region_3": {
        "region_3::oligo_2",
        "region_3::oligo_3",
        "region_3::oligo_4",
        "region_3::oligo_5",
    },
}

############################################
# Tests
############################################


def test_filter_exact_matches():
    # check that exact matches filters out a doubled sequence from the db
    oligo_database = OligoDatabase()
    exact_matches = ExactMatches()
    oligo_database.load_database(file_oligo_info_exact_matches)
    filtered_oligo_info_dict_match = exact_matches.apply(
        database=oligo_database.database,
        file_reference=None,
        n_jobs=n_jobs,
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
        coverage=coverage,
        word_size=word_size,
        percent_identity=percent_identity,
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
        coverage=coverage,
        word_size=word_size,
        percent_identity=percent_identity,
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


def test_cross_hybridization_filter_blast_bigger_region_policy(tmp_path):
    oligo_database = OligoDatabase()
    blastn = Blastn(
        tmp_path,
        word_size=word_size,
        percent_identity=percent_identity,
        strand="minus",
        coverage=coverage,
    )
    ch_filter = CrossHybridizationFilter(
        specificity_filter=blastn,
        policy=remove_by_bigger_region_policy,
        dir_cross_hybridization=None,
    )

    oligo_database.load_database(file_database=file_oligo_info_cross_hybridization)
    oligo_database_fasta = oligo_database.write_fasta_from_database(
        os.path.join(tmp_path, "oligo_DB_cross_hybridization")
    )

    filtered_database = ch_filter.apply(
        database=oligo_database.database,
        file_reference=oligo_database_fasta,
        n_jobs=None,
    )

    filtered_oligos = {
        key: {key_2 for key_2 in list(filtered_database[key].keys())}
        for key in list(filtered_database.keys())
    }
    assert (
        expected_oligos_bigger_region == filtered_oligos
    ), f"The cross-hybridization filter didn't return the expected oligos. \n\nExpected:\n{expected_oligos_bigger_region}\n\nGot:\n{filtered_oligos}"


def test_cross_hybridization_filter_bowtie_bigger_region_policy(tmp_path):
    oligo_database = OligoDatabase()
    bowtie_filter = Bowtie(tmp_path, num_mismatches, mismatch_region, strand="minus")
    ch_filter = CrossHybridizationFilter(
        specificity_filter=bowtie_filter,
        policy=remove_by_bigger_region_policy,
        dir_cross_hybridization=None,
    )

    oligo_database.load_database(file_database=file_oligo_info_cross_hybridization)
    oligo_database_fasta = oligo_database.write_fasta_from_database(
        os.path.join(tmp_path, "oligo_DB_cross_hybridization")
    )

    filtered_database = ch_filter.apply(
        database=oligo_database.database,
        file_reference=oligo_database_fasta,
        n_jobs=None,
    )

    filtered_oligos = {
        key: {key_2 for key_2 in list(filtered_database[key].keys())}
        for key in list(filtered_database.keys())
    }
    assert (
        expected_oligos_bigger_region == filtered_oligos
    ), f"The cross-hybridization filter didn't return the expected oligos. \n\nExpected:\n{expected_oligos_bigger_region}\n\nGot:\n{filtered_oligos}"


def test_cross_hybridization_filter_blast_degree_policy(tmp_path):
    oligo_database = OligoDatabase()
    blastn = Blastn(
        tmp_path,
        word_size=word_size,
        percent_identity=percent_identity,
        strand="minus",
        coverage=coverage,
    )
    ch_filter = CrossHybridizationFilter(
        specificity_filter=blastn,
        policy=remove_by_degree_policy,
        dir_cross_hybridization=None,
    )

    oligo_database.load_database(file_database=file_oligo_info_cross_hybridization)
    oligo_database_fasta = oligo_database.write_fasta_from_database(
        os.path.join(tmp_path, "oligo_DB_cross_hybridization")
    )

    filtered_database = ch_filter.apply(
        database=oligo_database.database,
        file_reference=oligo_database_fasta,
        n_jobs=None,
    )

    filtered_oligos = {
        key: {key_2 for key_2 in list(filtered_database[key].keys())}
        for key in list(filtered_database.keys())
    }
    assert (
        expected_oligos_degree == filtered_oligos
    ), f"The cross-hybridization filter didn't return the expected oligos. \n\nExpected:\n{expected_oligos_degree}\n\nGot:\n{filtered_oligos}"


def test_cross_hybridization_filter_bowtie_degree_policy(tmp_path):
    oligo_database = OligoDatabase()
    bowtie_filter = Bowtie(tmp_path, num_mismatches, mismatch_region, strand="minus")
    ch_filter = CrossHybridizationFilter(
        specificity_filter=bowtie_filter,
        policy=remove_by_degree_policy,
        dir_cross_hybridization=None,
    )

    oligo_database.load_database(file_database=file_oligo_info_cross_hybridization)
    oligo_database_fasta = oligo_database.write_fasta_from_database(
        os.path.join(tmp_path, "oligo_DB_cross_hybridization")
    )

    filtered_database = ch_filter.apply(
        database=oligo_database.database,
        file_reference=oligo_database_fasta,
        n_jobs=None,
    )

    filtered_oligos = {
        key: {key_2 for key_2 in list(filtered_database[key].keys())}
        for key in list(filtered_database.keys())
    }
    assert (
        expected_oligos_degree == filtered_oligos
    ), f"The cross-hybridization filter didn't return the expected oligos. \n\nExpected:\n{expected_oligos_degree}\n\nGot:\n{filtered_oligos}"
