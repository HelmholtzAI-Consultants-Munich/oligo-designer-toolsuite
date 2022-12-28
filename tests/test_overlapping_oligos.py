import shutil

import pandas as pd
import pytest

from oligo_designer_toolsuite.IO import CustomOligoDB
from oligo_designer_toolsuite.oligo_efficiency import (
    PadlockProbeScoring,
    PadlockSetScoring,
)
from oligo_designer_toolsuite.oligo_selection import (
    ProbesetGenerator,
    padlock_heuristic_selection,
)


@pytest.fixture()
def oligos_database():
    # database = NcbiOligoDB(probe_length_min=30, probe_length_max=40, filters=cls.filters, dir_output='tests/output')

    # If the anotation and fasta files are already saved on the machine, it is possible to direclty use them
    # instead of downloading them again.
    dir_annotation = "/home/francesco/Desktop/Work/NCBI"
    annotation = dir_annotation + "/GCF_000001405.40_GRCh38.p14_genomic.gtf"
    sequence = dir_annotation + "/GCF_000001405.40_GRCh38.p14_genomic.fna"
    database = CustomOligoDB(
        probe_length_min=30,
        probe_length_max=40,
        species="unknown",
        genome_assembly="unknown",
        annotation_release="unknown",
        annotation_source="unknown",
        file_annotation=annotation,
        file_sequence=sequence,
        dir_output="tests/output",
    )
    database.read_oligos_DB(
        format="tsv", file_oligos_DB_tsv="tests/data/probes_info.tsv"
    )

    yield database

    shutil.rmtree("tests/output")
    del database


@pytest.fixture()
def probeset_generator():
    padlock_scoring = PadlockProbeScoring(
        Tm_min=52,
        Tm_opt=60,
        Tm_max=67,
        GC_content_min=40,
        GC_content_opt=50,
        GC_content_max=60,
    )
    set_scoring = PadlockSetScoring()
    probeset_generator = ProbesetGenerator(
        probeset_size=5,
        min_probeset_size=2,
        probes_scoring=padlock_scoring,
        set_scoring=set_scoring,
        heurustic_selection=padlock_heuristic_selection,
    )
    return probeset_generator


# check we obtain the same result
def test_probesets_generation(probeset_generator, oligos_database):
    oligos_database = probeset_generator.apply(database=oligos_database, n_sets=100)
    for gene in oligos_database.probesets.keys():
        computed_sets = oligos_database.probesets[gene]
        computed_sets.drop(columns=["probeset_id"], inplace=True)
        true_sets = pd.read_csv(
            f"tests/data/probesets/ranked_probesets_{gene}.txt",
            sep="\t",
            index_col=0,
            float_precision="round_trip",
        )
        # additional sorting to guarantee the order is the same
        true_sets.sort_values(by=list(true_sets.columns), inplace=True)
        true_sets.reset_index(inplace=True, drop=True)
        computed_sets.sort_values(by=list(computed_sets.columns), inplace=True)
        computed_sets.reset_index(inplace=True, drop=True)
        print(true_sets)
        print(computed_sets)
        assert true_sets.equals(
            computed_sets
        ), f"Sets for {gene} are not computed correctly!"


# check the overlapping matrix is created correctly with 2 probes given in input
def test_nonoverlapping_matrix_ovelapping_probes(probeset_generator):
    # generate sinthetic probes that overlap
    probes = {
        "A_0": {"start": [10, 50], "end": [15, 55]},
        "A_1": {"start": [20, 53], "end": [25, 58]},
    }
    computed_matrix = probeset_generator._get_overlapping_matrix(probes=probes)
    true_matrix = pd.DataFrame(
        data=[[0, 0], [0, 0]], columns=["A_0", "A_1"], index=["A_0", "A_1"]
    )
    assert true_matrix.equals(
        computed_matrix
    ), "overlapping matrix for two overlapping probes wrongly computed"


def test_nonoverlapping_matrix_for_nonovelapping_probes(probeset_generator):
    # generate sinthetic probes that overlap
    probes = {
        "A_0": {"start": [10, 50], "end": [15, 55]},
        "A_1": {"start": [20, 35], "end": [25, 40]},
    }
    computed_matrix = probeset_generator._get_overlapping_matrix(probes=probes)
    true_matrix = pd.DataFrame(
        data=[[0, 1], [1, 0]], columns=["A_0", "A_1"], index=["A_0", "A_1"]
    )
    assert true_matrix.equals(
        computed_matrix
    ), "overlapping matrix for two non-overlapping probes wrongly computed"


def test_non_overlapping_sets(probeset_generator):
    probes = {
        "A_0": {"melting_temperature": 60, "GC_content": 50},
        "A_1": {"melting_temperature": 60, "GC_content": 55},
        "A_2": {"melting_temperature": 67, "GC_content": 55},
        "A_3": {"melting_temperature": 67, "GC_content": 60},
        "A_4": {"melting_temperature": 60, "GC_content": 60},
    }
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
            "probeset_id",
            "probe_0",
            "probe_1",
            "probe_2",
            "set_score_0",
            "set_score_1",
        ],
    )
    _, computed_sets, _ = probeset_generator._get_non_overlapping_sets(
        probes, overlapping_matrix, 100
    )
    assert true_sets.equals(computed_sets), "Sets are not computed correctly"
