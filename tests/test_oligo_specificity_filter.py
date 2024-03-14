############################################
# imports
############################################

import pytest

from oligo_designer_toolsuite.database import OligoDatabase, ReferenceDatabase
from oligo_designer_toolsuite.oligo_specificity_filter import (
    BlastNFilter,
    BlastNSeedregionLigationsiteFilter,
    Bowtie2Filter,
    BowtieFilter,
    CrossHybridizationFilter,
    ExactMatchFilter,
    RemoveByDegreePolicy,
    RemoveByLargerRegionPolicy,
)

############################################
# Global Parameters
############################################

# Files
file_database_oligos_exact_match = "data/tests/databases/database_oligos_exactmatch.tsv"
file_database_oligos_match = "data/tests/databases/database_oligos_match.tsv"
file_database_oligos_nomatch = "data/tests/databases/database_oligos_nomatch.tsv"

file_database_oligos_ligation_match = "data/tests/databases/database_oligos_ligation_match.tsv"
file_database_oligos_ligation_nomatch = "data/tests/databases/database_oligos_ligation_nomatch.tsv"

file_database_oligos_crosshyb = "data/tests/databases/database_oligos_crosshybridization.tsv"

file_database_reference = "data/tests/databases/database_reference.fna"
file_database_reference_ligation = "data/tests/databases/database_reference_ligation.fna"

file_annotation_ncbi = "data/annotations/custom_GCF_000001405.40_GRCh38.p14_genomic_chr16.gtf"
file_sequence_ncbi = "data/annotations/custom_GCF_000001405.40_GRCh38.p14_genomic_chr16.fna"

# Metadata
metadata_ncbi = {
    "files_source": "NCBI",
    "species": "Homo_sapiens",
    "annotation_release": "110",
    "genome_assembly": "GRCh38",
}

region_ids = ["AARS1", "DECR2", "FAM234A", "RHBDF1", "WASIR2"]

# Blast parameters
blast_search_parameters = {"perc_identity": 80, "strand": "plus", "word_size": 10}
blast_search_parameters_crosshyb = {"perc_identity": 80, "strand": "minus", "word_size": 10}

blast_hit_parameters = {"coverage": 50}
blast_hit_parameters_crosshyb = {"coverage": 50}

# Bowtie parameters
bowtie_search_parameters = {"-n": 3, "-l": 5}
bowtie_search_parameters_crosshyb = {"-n": 3, "-l": 5, "--nofw": ""}

bowtie2_search_parameters = {"-N": 0}

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
# Preprocessing
############################################


# tmp_path generates a temporary directory unique to the test invocation, created in the base temporary directory
# for more information see here: https://docs.pytest.org/en/6.2.x/tmpdir.html
@pytest.fixture(scope="session")
def oligo_database_exact_match(tmpdir_factory):
    base_temp = tmpdir_factory.getbasetemp()

    oligo_database_exact_match = OligoDatabase(
        min_oligos_per_region=2, write_regions_with_insufficient_oligos=True, dir_output=base_temp
    )
    oligo_database_exact_match.load_database(file_database_oligos_exact_match)

    return oligo_database_exact_match


@pytest.fixture(scope="session")
def oligo_database_match(tmpdir_factory):
    base_temp = tmpdir_factory.getbasetemp()

    oligo_database_match = OligoDatabase(
        min_oligos_per_region=2, write_regions_with_insufficient_oligos=True, dir_output=base_temp
    )
    oligo_database_match.load_database(file_database_oligos_match)

    return oligo_database_match


@pytest.fixture(scope="session")
def oligo_database_nomatch(tmpdir_factory):
    base_temp = tmpdir_factory.getbasetemp()

    oligo_database_nomatch = OligoDatabase(
        min_oligos_per_region=2, write_regions_with_insufficient_oligos=True, dir_output=base_temp
    )
    oligo_database_nomatch.load_database(file_database_oligos_nomatch)

    return oligo_database_nomatch


@pytest.fixture(scope="session")
def oligo_database_ligation_match(tmpdir_factory):
    base_temp = tmpdir_factory.getbasetemp()

    oligo_database_ligation_match = OligoDatabase(
        min_oligos_per_region=2, write_regions_with_insufficient_oligos=True, dir_output=base_temp
    )
    oligo_database_ligation_match.load_database(file_database_oligos_ligation_match)

    return oligo_database_ligation_match


@pytest.fixture(scope="session")
def oligo_database_ligation_nomatch(tmpdir_factory):
    base_temp = tmpdir_factory.getbasetemp()

    oligo_database_ligation_nomatch = OligoDatabase(
        min_oligos_per_region=2, write_regions_with_insufficient_oligos=True, dir_output=base_temp
    )
    oligo_database_ligation_nomatch.load_database(file_database_oligos_ligation_nomatch)

    return oligo_database_ligation_nomatch


@pytest.fixture(scope="session")
def oligo_database_crosshyb(tmpdir_factory):
    base_temp = tmpdir_factory.getbasetemp()

    oligo_database_crosshyb = OligoDatabase(
        min_oligos_per_region=2, write_regions_with_insufficient_oligos=True, dir_output=base_temp
    )
    oligo_database_crosshyb.load_database(file_database_oligos_crosshyb)

    return oligo_database_crosshyb


@pytest.fixture(scope="session")
def reference_database(tmpdir_factory):
    base_temp = tmpdir_factory.getbasetemp()

    reference_database = ReferenceDatabase(dir_output=base_temp)
    reference_database.load_sequences_fom_fasta(file_fasta=file_database_reference, database_overwrite=True)

    return reference_database


@pytest.fixture(scope="session")
def reference_database_ligation(tmpdir_factory):
    base_temp = tmpdir_factory.getbasetemp()

    reference_database_ligation = ReferenceDatabase(dir_output=base_temp)
    reference_database_ligation.load_sequences_fom_fasta(
        file_fasta=file_database_reference_ligation, database_overwrite=True
    )

    return reference_database_ligation


############################################
# Testing
############################################


def test_exact_match_filter(tmp_path, oligo_database_exact_match):
    sequence_type = "oligo"

    exactmatch_filter = ExactMatchFilter()
    res = exactmatch_filter.apply(sequence_type, oligo_database_exact_match, 2)

    assert (
        "WASH7P::2" not in res.database["WASH7P"].keys()
    ), "A matching oligo has not been filtered from exact matches!"
    assert (
        "AGRN::1" not in res.database["AGRN"].keys()
    ), "A matching oligo has not been filtered from exact mathces!"


def test_blastn_filter_match(tmp_path, oligo_database_match, reference_database):
    sequence_type = "target"

    blast_filter = BlastNFilter(blast_search_parameters, blast_hit_parameters, dir_output=tmp_path)
    res = blast_filter.apply(sequence_type, oligo_database_match, 2, reference_database)

    assert (
        "WASH7P::1" not in res.database["WASH7P"].keys()
    ), "A matching oligo has not been filtered by Blast!"


def test_blastn_filter_nomatch(tmp_path, oligo_database_nomatch, reference_database):
    sequence_type = "target"

    blast_filter = BlastNFilter(blast_search_parameters, blast_hit_parameters, dir_output=tmp_path)
    res = blast_filter.apply(sequence_type, oligo_database_nomatch, 2, reference_database)

    assert "AGRN::1" in res.database["AGRN"].keys(), "A non matching oligo has been filtered by Blast!"


def test_bowtie_filter_match(tmp_path, oligo_database_match, reference_database):
    sequence_type = "target"

    bowtie_filter = BowtieFilter(bowtie_search_parameters, dir_output=tmp_path)
    res = bowtie_filter.apply(sequence_type, oligo_database_match, 2, reference_database)

    assert (
        "WASH7P::1" not in res.database["WASH7P"].keys()
    ), "A matching oligo has not been filtered by Bowtie!"


def test_bowtie_filter_nomatch(tmp_path, oligo_database_nomatch, reference_database):
    sequence_type = "target"

    bowtie_filter = BowtieFilter(bowtie_search_parameters, dir_output=tmp_path)
    res = bowtie_filter.apply(sequence_type, oligo_database_nomatch, 2, reference_database)

    assert "AGRN::1" in res.database["AGRN"].keys(), "A non matching oligo has been filtered by Bowtie!"


def test_bowtie2_filter_match(tmp_path, oligo_database_match, reference_database):
    sequence_type = "target"

    bowtie_filter = Bowtie2Filter(bowtie2_search_parameters, dir_output=tmp_path)
    res = bowtie_filter.apply(sequence_type, oligo_database_match, 2, reference_database)

    assert (
        "WASH7P::1" not in res.database["WASH7P"].keys()
    ), "A matching oligo has not been filtered by Bowtie2!"


def test_bowtie2_filter_nomatch(tmp_path, oligo_database_nomatch, reference_database):
    sequence_type = "target"

    bowtie_filter = Bowtie2Filter(bowtie2_search_parameters, dir_output=tmp_path)
    res = bowtie_filter.apply(sequence_type, oligo_database_nomatch, 2, reference_database)

    assert "AGRN::1" in res.database["AGRN"].keys(), "A non matching oligo has been filtered by Bowtie2!"


def test_blastn_ligation_filter_match(tmp_path, oligo_database_ligation_match, reference_database_ligation):
    sequence_type = "target"
    seedregion_size = 10

    blast_ligation_filter = BlastNSeedregionLigationsiteFilter(
        seedregion_size, blast_search_parameters, blast_hit_parameters, dir_output=tmp_path
    )
    res = blast_ligation_filter.apply(
        sequence_type, oligo_database_ligation_match, 2, reference_database_ligation
    )

    assert (
        "WASH7P::1" not in res.database["WASH7P"].keys()
    ), "A matching oligo has not been filtered by Blast!"


def test_blastn_ligation_filter_nomatch(
    tmp_path, oligo_database_ligation_nomatch, reference_database_ligation
):
    sequence_type = "target"
    seedregion_size = 5

    blast_ligation_filter = BlastNSeedregionLigationsiteFilter(
        seedregion_size, blast_search_parameters, blast_hit_parameters, dir_output=tmp_path
    )
    res = blast_ligation_filter.apply(
        sequence_type, oligo_database_ligation_nomatch, 2, reference_database_ligation
    )

    assert "AGRN::1" in res.database["AGRN"].keys(), "A non matching oligo has been filtered by Blast!"


def test_crosshyb_filter_exactmatch_bigger_region_policy(tmp_path, oligo_database_exact_match):
    exactmatch_filter = ExactMatchFilter()
    policy = RemoveByLargerRegionPolicy()

    sequence_type = "oligo"

    cross_hyb_filter = CrossHybridizationFilter(policy, exactmatch_filter, tmp_path)
    res = cross_hyb_filter.apply(sequence_type, oligo_database_exact_match, 2)

    assert (
        "WASH7P::1" in res.database["WASH7P"].keys()
    ), "A non matching oligo has been filtered by exact matches!"
    assert (
        "WASH7P::2" not in res.database["WASH7P"].keys()
    ), "A matching oligo has not been filtered by exact mathces!"
    assert (
        "AGRN::1" in res.database["AGRN"].keys()
    ), "A non matching oligo has been filtered by exact matches!"


def test_crosshyb_filter_blast_bigger_region_policy(tmp_path, oligo_database_crosshyb):
    blast_filter = BlastNFilter(
        blast_search_parameters_crosshyb, blast_hit_parameters_crosshyb, dir_output=tmp_path
    )
    policy = RemoveByLargerRegionPolicy()

    sequence_type = "oligo"

    cross_hyb_filter = CrossHybridizationFilter(policy, blast_filter, tmp_path)
    res = cross_hyb_filter.apply(sequence_type, oligo_database_crosshyb, 2)

    filtered_oligos = {
        key: {key_2 for key_2 in list(res.database[key].keys())} for key in list(res.database.keys())
    }
    assert (
        expected_oligos_bigger_region == filtered_oligos
    ), f"The cross-hybridization filter didn't return the expected oligos. \n\nExpected:\n{expected_oligos_bigger_region}\n\nGot:\n{filtered_oligos}"


def test_crosshyb_filter_blast_degree_policy(tmp_path, oligo_database_crosshyb):
    blast_filter = BlastNFilter(
        blast_search_parameters_crosshyb, blast_hit_parameters_crosshyb, dir_output=tmp_path
    )
    policy = RemoveByDegreePolicy()

    sequence_type = "oligo"

    cross_hyb_filter = CrossHybridizationFilter(policy, blast_filter, tmp_path)
    res = cross_hyb_filter.apply(sequence_type, oligo_database_crosshyb, 2)

    filtered_oligos = {
        key: {key_2 for key_2 in list(res.database[key].keys())} for key in list(res.database.keys())
    }
    assert (
        expected_oligos_degree == filtered_oligos
    ), f"The cross-hybridization filter didn't return the expected oligos. \n\nExpected:\n{expected_oligos_degree}\n\nGot:\n{filtered_oligos}"


def test_crosshyb_filter_bowtie_bigger_region_policy(tmp_path, oligo_database_crosshyb):
    bowtie_filter = BowtieFilter(bowtie_search_parameters_crosshyb, dir_output=tmp_path)
    policy = RemoveByLargerRegionPolicy()

    sequence_type = "oligo"

    cross_hyb_filter = CrossHybridizationFilter(policy, bowtie_filter, tmp_path)
    res = cross_hyb_filter.apply(sequence_type, oligo_database_crosshyb, 2)

    filtered_oligos = {
        key: {key_2 for key_2 in list(res.database[key].keys())} for key in list(res.database.keys())
    }
    assert (
        expected_oligos_bigger_region == filtered_oligos
    ), f"The cross-hybridization filter didn't return the expected oligos. \n\nExpected:\n{expected_oligos_bigger_region}\n\nGot:\n{filtered_oligos}"


def test_crosshyb_filter_bowtie_degree_policy(tmp_path, oligo_database_crosshyb):
    bowtie_filter = BowtieFilter(bowtie_search_parameters_crosshyb, dir_output=tmp_path)
    policy = RemoveByDegreePolicy()

    sequence_type = "oligo"

    cross_hyb_filter = CrossHybridizationFilter(policy, bowtie_filter, tmp_path)
    res = cross_hyb_filter.apply(sequence_type, oligo_database_crosshyb, 2)

    filtered_oligos = {
        key: {key_2 for key_2 in list(res.database[key].keys())} for key in list(res.database.keys())
    }
    assert (
        expected_oligos_degree == filtered_oligos
    ), f"The cross-hybridization filter didn't return the expected oligos. \n\nExpected:\n{expected_oligos_degree}\n\nGot:\n{filtered_oligos}"
