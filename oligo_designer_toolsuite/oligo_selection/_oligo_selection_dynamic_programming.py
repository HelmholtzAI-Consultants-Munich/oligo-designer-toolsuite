############################################
# imports
############################################

import bisect
from math import ceil, inf

import numpy as np
import pandas as pd

from oligo_designer_toolsuite.database import OligoDatabase
from oligo_designer_toolsuite.oligo_efficiency_filter import (
    AverageSetScoring,
    LowestSetScoring,
    OligoScoring,
    SetScoringBase,
)

from ._oligo_selection_base import BaseOligoSelection

############################################
# Oligoset Selection Classes
############################################


class DynamicProgrammingOligoSelection(BaseOligoSelection):

    def __init__(
        self,
        oligos_scoring: OligoScoring,
        set_scoring: SetScoringBase,
        set_size_opt: int,
        set_size_min: int,
        distance_between_oligos: int,
        jaccard_opt: float,
    ) -> None:
        """Constructor for the IndependentSetsOligoSelection class."""
        self.oligos_scoring = oligos_scoring
        self.set_scoring = set_scoring
        self.set_size_opt = set_size_opt
        self.set_size_min = set_size_min
        self.distance_between_oligos = distance_between_oligos
        self.jaccard_opt = jaccard_opt

    def _get_oligo_sets_for_region(
        self,
        oligo_database: OligoDatabase,
        sequence_type: str,
        region_id: str,
        n_sets: int,
    ) -> None:
        """
        Generates oligo sets for a specific region by using dynamic programming.
        """

        database_region: dict[str, dict] = oligo_database.database[region_id]
        oligo_ids: list[str] = list(database_region.keys())
        n = len(oligo_ids)

        # Treat oligos with multiple positions as very long oligos, stretching from their first occurence to their last one
        # Multiple positions mostly due to oligos spanning exons (i.e. on exon-exon junctions)
        # These are treated as incompatible with oligos within the bridged introns
        start_points: list[int] = [
            min(min(starts) for starts in database_region[id]["start"]) for id in oligo_ids
        ]
        end_points: list[int] = [max(max(ends) for ends in database_region[id]["end"]) for id in oligo_ids]

        # Prepare reordering
        order = list(range(n))
        order.sort(key=lambda i: end_points[i])
        # Used later to map results to original oligos
        reordered_oligo_ids = [""] * n
        for new_i, i in enumerate(order):
            reordered_oligo_ids[new_i] = oligo_ids[i]

        # Sort intervals asc. by end points
        start_points = [start_points[i] for i in order]
        end_points = [end_points[i] for i in order]

        # Calculate scores
        scores_series = self.oligos_scoring.apply(
            oligo_database, region_id, reordered_oligo_ids, sequence_type
        )

        # Invert scores because the lower the better, but we maximize
        scores = [-s for s in scores_series]

        skip_indices: set[int] = set()  # set of indices already used in an oligoset
        oligosets = []
        for set_idx in range(n_sets):
            # Calculate optimal set for remaining oligos
            nodes = self._get_oligo_set_for_subset(start_points, end_points, scores, skip_indices)
            if not nodes:
                break

            skip_indices.update(self._select_nodes_to_skip(nodes, scores))
            oligosets.append({reordered_oligo_ids[i] for i in nodes})

        oligo_database.oligosets[region_id] = self._oligosets_to_dataframe(
            oligosets, scores_series, self.set_size_opt
        )

        # Remove oligos from database that are not part of oligosets
        oligos_keep = set().union(*oligosets)
        for oligo_id in oligo_ids:
            if oligo_id not in oligos_keep:
                del database_region[oligo_id]
        oligo_database.database[region_id] = database_region

    def _select_nodes_to_skip(self, oligoset: list[int], scores: list[int]) -> set[int]:
        """Chooses nodes that should be ignored for finding additional oligosets.

        Selects the worst k nodes, where `k = n_nodes * (1-self.jaccard_opt)`.
        This ensures that new oligosets can not have a Jaccard similarity greater
        than `self.jaccard_opt` to any previously chosen oligoset.
        """
        k_nodes_to_remove = max(ceil((1 - self.jaccard_opt) * len(oligoset)), 1)
        sorted_oligoset = sorted(oligoset, key=lambda i: scores[i])
        return set(sorted_oligoset[:k_nodes_to_remove])

    def _get_oligo_set_for_subset(
        self,
        start_points: list[int],
        end_points: list[int],
        scores: list[float],
        skip_indices: set[int],
    ) -> list[int]:
        total_n = len(start_points)

        remaining_indices = [i for i in range(total_n) if i not in skip_indices]
        start_points = [start_points[i] for i in remaining_indices]
        end_points = [end_points[i] for i in remaining_indices]
        scores = [scores[i] for i in remaining_indices]

        selected_subset: list[int]
        match self.set_scoring:
            case LowestSetScoring():
                selected_subset = self._maximum_lowest_weight_oligo_selection(
                    start_points, end_points, scores, self.set_size_min, self.set_size_opt
                )
            case AverageSetScoring():
                selected_subset = self._maximum_weight_sum_oligo_selection(
                    start_points, end_points, scores, self.set_size_min, self.set_size_opt
                )
            case _:
                raise NotImplementedError()

        # Map back to original indices
        #   e.g. if skip_indices == {0} the above function may include 0 in its result,
        #   meaning the first element of the remaining nodes (i.e. 1)
        return [remaining_indices[i] for i in selected_subset]

    def _maximum_weight_sum_oligo_selection(
        self,
        start_points: list[int],
        end_points: list[int],
        scores: list[float],
        minimum_set_size: int,
        maximum_set_size: int,
    ) -> list[int]:
        """Select subset such that the sum of its weights is maximum.

        For a fixed set size, this is equivalent to finding the subset of maximum average weight.
        Input intervals must be sorted by end points.
        """
        n = len(start_points)

        if (minimum_set_size > maximum_set_size) or (minimum_set_size > n) or (minimum_set_size < 1):
            return []

        def reconstruct_result(k: int) -> list[int]:
            # reconstruct nodes in O(klogn+n)
            nodes = []

            i = n - 1
            for j in reversed(range(1, k + 1)):
                bisect_idx = (
                    bisect.bisect(used_nodes[j], i) - 1
                )  # either the index of i or of the next-smallest in used_nodes[j]
                node = used_nodes[j][bisect_idx]
                nodes.append(node)
                i = preds[node]

            nodes.reverse()
            return nodes

        suboptimal = False
        used_nodes: list[list[int]] = [[] for _ in range(maximum_set_size + 1)]

        # preprocessing: find predecessors O(nlogn)
        preds = (np.searchsorted(end_points, start_points, side="right") - 1).tolist()

        # (sum_of_scores, lowest_score)
        dp_prev: list[tuple[float, float]] = list(zip(scores, scores))
        dp_cur: list[tuple[float, float]] = [(-inf, -inf)] * n

        # add one last element to serve as fallback for intervals that don't have predecessor
        # if there's no predecessor, preds[i] will be -1 and using that as an index will go to this last element
        # this lets us avoid an additional branch in the inner loop
        dp_prev.append((-inf, -inf))
        dp_cur.append((-inf, -inf))

        # initialization for j=1 O(n)
        used_nodes[1].append(0)
        for i in range(1, n):
            if (last := dp_prev[i - 1]) > dp_prev[i]:  # TODO: >= here for consistency with below
                dp_prev[i] = last
            else:
                used_nodes[1].append(i)

        # processing O(n*k_max)
        for j in range(2, maximum_set_size + 1):
            dp_cur[0] = (-inf, -inf)
            running_best = (-inf, -inf)
            for i in range(1, n):
                score = scores[i]

                pred_sum, pred_lowest = dp_prev[preds[i]]
                # if no predecessor exists or it doesn't have a value for j-1, it's (-inf, -inf)
                if score >= pred_lowest:
                    new_val: tuple[float, float] = (pred_sum + score, pred_lowest)
                else:
                    new_val = (pred_sum + score, score)

                if running_best >= new_val:
                    dp_cur[i] = running_best
                else:
                    dp_cur[i] = new_val
                    running_best = new_val
                    used_nodes[j].append(i)

            if dp_cur[-2][0] == -inf:
                suboptimal = True
                break
            dp_prev, dp_cur = dp_cur, dp_prev

        # j-1 is the last successful loop if it's suboptimal
        result_size = maximum_set_size if not suboptimal else j - 1
        # check whether minimum set size achieved
        if result_size < minimum_set_size:
            return []

        return reconstruct_result(result_size)

    def _maximum_lowest_weight_oligo_selection(
        self,
        start_points: list[int],
        end_points: list[int],
        scores: list[float],
        minimum_set_size: int,
        maximum_set_size: int,
    ) -> list[int]:
        """Select subset such that lowest weight is maximum.

        Input intervals must be sorted by end points.
        """
        n = len(start_points)

        if (minimum_set_size > maximum_set_size) or (minimum_set_size > n) or (minimum_set_size < 1):
            return []

        def reconstruct_result(k: int) -> list[int]:
            # reconstruct nodes in O(klogn+n)
            nodes = []

            i = n - 1
            for j in reversed(range(1, k + 1)):
                bisect_idx = (
                    bisect.bisect(used_nodes[j], i) - 1
                )  # either the index of i or of the next-smallest in used_nodes[j]
                node = used_nodes[j][bisect_idx]
                nodes.append(node)
                i = preds[node]

            nodes.reverse()
            return nodes

        suboptimal = False
        used_nodes: list[list[int]] = [[] for _ in range(maximum_set_size + 1)]

        # preprocessing: find predecessors O(nlogn)
        preds = (np.searchsorted(end_points, start_points, side="right") - 1).tolist()

        # (lowest_score, sum_of_scores)
        dp_prev: list[tuple[float, float]] = list(zip(scores, scores))
        dp_cur: list[tuple[float, float]] = [(-inf, -inf)] * n

        # add one last element to serve as fallback for intervals that don't have predecessor
        # if there's no predecessor, preds[i] will be -1 and using that as an index will go to this last element
        # this lets us avoid an additional branch in the inner loop
        dp_prev.append((-inf, -inf))
        dp_cur.append((-inf, -inf))

        # initialization for j=1 O(n)
        used_nodes[1].append(0)
        for i in range(1, n):
            if (last := dp_prev[i - 1]) > dp_prev[i]:
                dp_prev[i] = last
            else:
                used_nodes[1].append(i)

        # processing O(n*k_max)
        for j in range(2, maximum_set_size + 1):
            dp_cur[0] = (-inf, -inf)
            running_best = (-inf, -inf)
            for i in range(1, n):
                score = scores[i]

                pred_lowest, pred_sum = dp_prev[preds[i]]
                # if no predecessor exists or it doesn't have a value for j-1, it's (-inf, -inf)
                if score >= pred_lowest:
                    new_val: tuple[float, float] = (pred_lowest, pred_sum + score)
                else:
                    new_val = (score, pred_sum + score)

                if running_best >= new_val:
                    dp_cur[i] = running_best
                else:
                    dp_cur[i] = new_val
                    running_best = new_val
                    used_nodes[j].append(i)

            if dp_cur[-2][0] == -inf:
                suboptimal = True
                break
            dp_prev, dp_cur = dp_cur, dp_prev

        # j-1 is the last successful loop if it's suboptimal
        result_size = maximum_set_size if not suboptimal else j - 1
        # check whether minimum set size achieved
        if result_size < minimum_set_size:
            return []

        return reconstruct_result(result_size)

    def _oligosets_to_dataframe(
        self,
        oligosets: list[set[str]],
        scores: pd.Series,
        maximum_set_size: int,
    ) -> pd.DataFrame:
        """
        Convert selected oligo sets and their scores into a tidy `DataFrame`.

        Each row corresponds to one oligo set, with columns:
        - `oligoset_id`: an integer identifier
        - `oligo_0`, ..., `oligo_n`: oligo IDs in the set
        - one column per set score in `set_scoring.score_names`

        :param oligosets: Dictionary mapping oligo tuples to their set score dicts.
        :type oligosets: dict[tuple[str, ...], dict[str, float]]
        :return: DataFrame with one row per oligo set.
        :rtype: pd.DataFrame
        """
        if not oligosets:
            return pd.DataFrame()

        rows = []
        for oligoset in oligosets:
            oligoset_list = list(oligoset)
            set_scores = scores[oligoset_list].round(4)
            average_set_score = set_scores.mean()
            lowest_set_score = set_scores.max()
            sum_set_score = set_scores.sum()

            # Pad oligoset to maximum_set_size if len(oligoset) < maximum_set_size
            padded_oligo_ids = oligoset_list + [None] * (maximum_set_size - len(oligoset_list))
            rows.append([*padded_oligo_ids, average_set_score, lowest_set_score, sum_set_score])

        columns = [f"oligo_{i}" for i in range(maximum_set_size)] + [
            "set_score_average",
            "set_score_worst",
            "set_score_sum",
        ]

        df = pd.DataFrame(rows, columns=columns)
        df.insert(0, "oligoset_id", range(len(df)))
        return df
