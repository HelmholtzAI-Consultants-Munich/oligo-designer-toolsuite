############################################
# imports
############################################
from pandas import Series

from oligo_designer_toolsuite.oligo_efficiency_filter import (
    AverageSetScoring,
    MaxSetScoring,
    PadlockOligoScoring,
    PadlockSetScoring,
)

############################################
# Global Parameters
############################################

set = Series(data=[0, 1, 8, 5, 2, 6, 7, 3])
n = 5

############################################
# Tests
############################################


def test_padlock():
    score = PadlockSetScoring()
    oligoset = score.apply(set, n)
    assert oligoset == [0, 1, 4, 7, 3, 5, 11], "padlock scoring failed"


def test_avg():
    score = AverageSetScoring()
    oligoset = score.apply(set, n)
    assert oligoset == [0, 1, 4, 7, 3, 2.2], "Average scoring failed"


def test_max():
    score = MaxSetScoring()
    oligoset = score.apply(set, n)
    assert oligoset == [0, 1, 4, 7, 3, 5], "Max scoring failed"


def test_padlock_oligo():
    score = PadlockOligoScoring(0, 1, 2, 0.25, 0.5, 0.75)
    oligo = {"TmNN": 0.5, "GC_content": 0.55}
    oligo_score = score.scoring_function(oligo)
    assert abs(oligo_score - 0.7) < 1e-5, "Oligo padlock score failed!"
