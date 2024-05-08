############################################
# imports
############################################

import pandas as pd
import pytest

from oligo_designer_toolsuite.database import OligoDatabase
from oligo_designer_toolsuite.oligo_efficiency_filter import (
    PadlockOligoScoring,
    PadlockSetScoring,
)
from oligo_designer_toolsuite.oligo_selection import (
    OligosetGenerator,
    padlock_heuristic_selection,
)

############################################
# Global Parameters
############################################

annotation_file_ncbi = "data/annotations/custom_GCF_000001405.40_GRCh38.p14_genomic_chr16.gtf"
sequence_file_ncbi = "data/annotations/custom_GCF_000001405.40_GRCh38.p14_genomic_chr16.fna"

############################################
# Tests
############################################


@pytest.fixture(scope="session")
def oligos_database(tmpdir_factory):
    base_temp = tmpdir_factory.getbasetemp()
    metadata = {"species": "human", "annotation_release": "110"}
    database = OligoDatabase(
        dir_output=base_temp,
    )
    database.load_metadata(metadata=metadata)
    database.load_database("data/tests/oligo_selection/oligos_info.tsv")

    yield database


@pytest.fixture()
def oligoset_generator():
    padlock_scoring = PadlockOligoScoring(
        Tm_min=52,
        Tm_opt=60,
        Tm_max=67,
        GC_content_min=40,
        GC_content_opt=50,
        GC_content_max=60,
    )
    set_scoring = PadlockSetScoring()
    oligoset_generator = OligosetGenerator(
        oligoset_size=5,
        min_oligoset_size=2,
        oligos_scoring=padlock_scoring,
        set_scoring=set_scoring,
        heurustic_selection=padlock_heuristic_selection,
    )
    return oligoset_generator


@pytest.fixture()
def padlock_scoring():
    padlock_scoring = PadlockOligoScoring(
        Tm_min=52,
        Tm_opt=60,
        Tm_max=67,
        GC_content_min=40,
        GC_content_opt=50,
        GC_content_max=60,
    )
    return padlock_scoring


# check we obtain the same result
def test_oligosets_generation(oligoset_generator, oligos_database):
    oligos_database = oligoset_generator.apply(oligo_database=oligos_database, n_sets=100, n_jobs=1)
    for gene in oligos_database.oligosets.keys():
        computed_sets = oligos_database.oligosets[gene]
        computed_sets.drop(columns=["oligoset_id"], inplace=True)
        true_sets = pd.read_csv(
            f"data/tests/oligo_selection/ranked_oligosets_{gene}.txt",
            sep="\t",
            index_col=0,
            float_precision="round_trip",
        )
        # additional sorting to guarantee the order is the same
        true_sets.sort_values(by=list(true_sets.columns), inplace=True)
        true_sets.reset_index(inplace=True, drop=True)
        computed_sets.sort_values(by=list(computed_sets.columns), inplace=True)
        computed_sets.reset_index(inplace=True, drop=True)
        pd.set_option("display.max_columns", 500)
        assert true_sets.equals(computed_sets), f"Sets for {gene} are not computed correctly!"


# check the overlapping matrix is created correctly with 2 oligos given in input
def test_nonoverlapping_matrix_ovelapping_oligos(oligoset_generator):
    # generate sinthetic oligos that overlap
    oligos = {
        "A_0": {"start": [10, 50], "end": [15, 55]},
        "A_1": {"start": [20, 53], "end": [25, 58]},
    }
    computed_matrix = oligoset_generator._get_overlapping_matrix(database_region=oligos)
    true_matrix = pd.DataFrame(data=[[0, 0], [0, 0]], columns=["A_0", "A_1"], index=["A_0", "A_1"])
    assert true_matrix.equals(
        computed_matrix
    ), "overlapping matrix for two overlapping oligos wrongly computed"


def test_nonoverlapping_matrix_for_nonovelapping_oligos(oligoset_generator):
    # generate sinthetic oligos that overlap
    oligos = {
        "A_0": {"start": [10, 50], "end": [15, 55]},
        "A_1": {"start": [20, 35], "end": [25, 40]},
    }
    computed_matrix = oligoset_generator._get_overlapping_matrix(database_region=oligos)
    true_matrix = pd.DataFrame(data=[[0, 1], [1, 0]], columns=["A_0", "A_1"], index=["A_0", "A_1"])
    assert true_matrix.equals(
        computed_matrix
    ), "overlapping matrix for two non-overlapping oligos wrongly computed"


def test_non_overlapping_sets(oligoset_generator, padlock_scoring):
    oligos = {
        "A_0": {"melting_temperature": 60, "GC_content": 50},
        "A_1": {"melting_temperature": 60, "GC_content": 55},
        "A_2": {"melting_temperature": 67, "GC_content": 55},
        "A_3": {"melting_temperature": 67, "GC_content": 60},
        "A_4": {"melting_temperature": 60, "GC_content": 60},
    }
    oligos, oligo_scores = padlock_scoring.apply(oligos)
    index = ["A_0", "A_1", "A_2", "A_3", "A_4"]
    data = [
        [0, 1, 1, 0, 0],
        [1, 0, 1, 0, 0],
        [1, 1, 0, 1, 1],
        [0, 0, 1, 0, 1],
        [0, 0, 1, 1, 0],
    ]
    overlapping_matrix = pd.DataFrame(data=data, index=index, columns=index)
    data = [[0, "A_0", "A_1", "A_2", 1.5, 2.0], [1, "A_4", "A_2", "A_3", 2.0, 4.5]]
    true_sets = pd.DataFrame(
        data=data,
        columns=[
            "oligoset_id",
            "oligo_0",
            "oligo_1",
            "oligo_2",
            "set_score_0",
            "set_score_1",
        ],
    )
    _, computed_sets, _ = oligoset_generator._get_non_overlapping_sets(
        oligos, overlapping_matrix, oligo_scores, 100
    )
    assert true_sets.equals(computed_sets), "Sets are not computed correctly"
