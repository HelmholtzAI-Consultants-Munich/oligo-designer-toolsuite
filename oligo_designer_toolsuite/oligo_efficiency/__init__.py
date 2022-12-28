"""
This module assignes a score to the probes and the probesets for their (on-target) efficiency.
"""


from ._probe_scoring import PadlockProbeScoring, ProbeScoringBase
from ._set_scoring import (
    AverageSetScoring,
    MaxSetScoring,
    PadlockSetScoring,
    SetScoringBase,
)

__all__ = [
    "ProbeScoringBase",
    "PadlockProbeScoring",
    "SetScoringBase",
    "PadlockSetScoring",
    "AverageSetScoring",
    "MaxSetScoring",
]
