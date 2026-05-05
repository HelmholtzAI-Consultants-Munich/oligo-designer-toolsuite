############################################
# imports
############################################

from typing import Any

from scipy.sparse import csr_matrix

from oligo_designer_toolsuite.database import OligoDatabase
from oligo_designer_toolsuite.oligo_efficiency_filter import BaseScorer

############################################
# Set Property Scorer Classes
############################################


class UniformDistanceScorer(BaseScorer):
    """
    Scores an oligo based on how uniformly spaced it is relative to other oligos in a candidate set.

    Uses the non-overlap matrix (pairwise distances) to compute an optimal spacing from the
    maximum distance and average oligo length, then penalizes deviation of the oligo's minimum
    distance to the set from that optimum. Lower scores indicate more uniform spacing.

    :param score_weight: Weight applied to the uniform-distance score.
    :type score_weight: float
    """

    def __init__(self, average_oligo_length: float, score_weight: float) -> None:
        """Constructor for the UniformDistanceScorer class."""
        self.average_oligo_length = average_oligo_length
        self.score_weight = score_weight

    def apply(  # type: ignore[override]
        self,
        oligo_database: OligoDatabase,
        region_id: str,
        oligo_id: str,
        sequence_type: str,
        *,
        non_overlap_matrix: csr_matrix,
        non_overlap_matrix_ids: list[str],
        set_oligo_ids: list[str],
        oligoset_size: int,
        **_kwargs: Any,
    ) -> float:
        """
        Compute the uniform-distance score for the given oligo with respect to the candidate set.

        :param oligo_database: The OligoDatabase instance containing oligonucleotide sequences and their associated properties.
        :type oligo_database: OligoDatabase
        :param region_id: Region ID to process.
        :type region_id: str
        :param oligo_id: The ID of the oligo for which the score is computed.
        :type oligo_id: str
        :param sequence_type: Type of sequence being processed.
        :type sequence_type: str
        :param non_overlap_matrix: Sparse matrix of pairwise distances (non-overlap) between oligos.
        :type non_overlap_matrix: csr_matrix
        :param non_overlap_matrix_ids: List of oligo IDs corresponding to the matrix indices.
        :type non_overlap_matrix_ids: list[str]
        :param set_oligo_ids: List of oligo IDs in the current candidate set.
        :type set_oligo_ids: list[str]
        :param oligoset_size: Target size of the oligo set.
        :type oligoset_size: int
        :return: Weighted score reflecting deviation from uniform spacing (lower is better).
        :rtype: float
        """
        if not set_oligo_ids or oligoset_size < 2 or self.score_weight == 0:
            return 0.0

        d_max = float(non_overlap_matrix.data.max())
        dist_opt = max(0.0, (d_max - (oligoset_size - 2) * self.average_oligo_length)) / (oligoset_size - 1)

        oligo_idx = non_overlap_matrix_ids.index(oligo_id)
        set_idxs = [non_overlap_matrix_ids.index(id) for id in set_oligo_ids]

        row = non_overlap_matrix[oligo_idx, set_idxs]
        if row.data.size == 0 or dist_opt == 0:
            d_min = 0.0
        else:
            d_min = float(row.data.min())
        score = abs(dist_opt - d_min) / dist_opt
        return score * self.score_weight
