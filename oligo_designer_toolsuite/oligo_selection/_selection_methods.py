############################################
# imports
############################################

from abc import abstractmethod
from typing import List

import networkx as nx
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

from oligo_designer_toolsuite.oligo_efficiency_filter import SetScoringBase

############################################
# Heuristic Selection Methods
############################################


class OligoSelectionPolicy:

    def __init__(
        self,
        set_size_opt: int,
        set_size_min: int,
        n_sets: int,
        ascending: bool,
        set_scoring: SetScoringBase,
        pre_filter: bool,
    ) -> None:
        self.opt_oligoset_size = set_size_opt
        self.min_oligoset_size = set_size_min
        self.n_sets = n_sets
        self.ascending = ascending
        self.set_scoring = set_scoring
        self.pre_filter = pre_filter

    def apply(
        self, oligos_scores: pd.Series, overlapping_matrix: csr_matrix, n_attempts: int
    ) -> pd.DataFrame:

        if self.pre_filter:
            oligos_scores, overlapping_matrix = self.pre_filter_oligos(oligos_scores, overlapping_matrix)

        if len(oligos_scores) > self.min_oligoset_size:
            oligosets = self.run_selection(oligos_scores, overlapping_matrix, n_attempts)
        else:
            oligosets = []

        # not clean but it works
        score_names = list(self.set_scoring.apply(pd.Series([], dtype=object), 0)[1].keys())

        return self._format_as_df(oligosets, score_names)

    def _format_as_df(self, oligosets: list, score_names: list) -> pd.DataFrame:

        if len(oligosets) > 0:
            oligosets_columns = [f"oligo_{i}" for i in range(self.opt_oligoset_size)] + score_names
            oligosets = pd.DataFrame(
                columns=oligosets_columns,
                data=oligosets,
            )
            # Sort oligosets by score
            oligosets.drop_duplicates(inplace=True, subset=oligosets.columns[:-1])
            oligosets.sort_values(
                list(oligosets.columns[self.opt_oligoset_size :]), ascending=self.ascending, inplace=True
            )
            oligosets = oligosets.head(self.n_sets)
            oligosets.reset_index(drop=True, inplace=True)
            oligosets.insert(0, "oligoset_id", oligosets.index)
            return oligosets
        else:
            return None

    def pre_filter_oligos(
        self, oligos_scores: pd.Series, overlapping_matrix: csr_matrix
    ) -> tuple[pd.Series, csr_matrix]:

        G = nx.from_scipy_sparse_array(overlapping_matrix)
        max_clique = nx.approximation.max_clique(G)
        nodes = []
        while len(max_clique) > self.min_oligoset_size:
            nodes += list(max_clique)
            G.remove_nodes_from(max_clique)
            max_clique = nx.approximation.max_clique(G)

        oligos_scores = oligos_scores.iloc[nodes]
        overlapping_matrix = overlapping_matrix[nodes, :][:, nodes]
        return oligos_scores, overlapping_matrix

    @abstractmethod
    def run_selection(
        self, oligos_scores: pd.Series, overlapping_matrix: csr_matrix, n_attempts: int
    ) -> list:

        pass


class GreedySelectionPolicy(OligoSelectionPolicy):
    def __init__(
        self,
        set_size_opt: int,
        set_size_min: int,
        n_sets: int,
        ascending: bool,
        set_scoring: SetScoringBase,
        score_criteria: str,
        penalty: float = 0.05,
    ) -> None:

        super().__init__(set_size_opt, set_size_min, n_sets, ascending, set_scoring)
        self.score_criteria = score_criteria
        self.penalty = penalty

    def run_selection(
        self, oligos_scores: pd.Series, overlapping_matrix: csr_matrix, n_attempts: int
    ) -> List[List[int]]:

        oligo_to_idx = {oligo: idx for idx, oligo in enumerate(oligos_scores.index)}
        selected_sets = set()
        adjusted_scores = oligos_scores.copy()

        for set_size in range(self.opt_oligoset_size, self.min_oligoset_size - 1, -1):
            attempts = 0
            while len(selected_sets) < self.n_sets and attempts < n_attempts:
                attempts += 1
                best_set = []
                available_oligos = adjusted_scores.sort_values(ascending=self.ascending).index.to_list()
                for _ in range(set_size):
                    best_oligo = None
                    best_score = None

                    for oligo in available_oligos:
                        if not all(
                            overlapping_matrix[oligo_to_idx[oligo], oligo_to_idx[best_set_oligo]]
                            for best_set_oligo in best_set
                        ):
                            continue

                        current_set = best_set + [oligo]
                        current_scores = adjusted_scores.loc[current_set]

                        _, scores = self.set_scoring.apply(current_scores, set_size)

                        score = scores[self.score_criteria]
                        better_than_best = (
                            best_score is None
                            or (self.ascending and score < best_score)
                            or (not self.ascending and score > best_score)
                        )
                        if better_than_best:
                            best_oligo = oligo
                            best_score = score

                            oligo_scores_report = oligos_scores[current_set]
                            _, best_report = self.set_scoring.apply(oligo_scores_report, set_size)

                    if best_oligo is not None:
                        best_set.append(best_oligo)
                        available_oligos.remove(best_oligo)
                        adjusted_scores[best_oligo] *= (
                            (1 + self.penalty) if self.ascending else (1 - self.penalty)
                        )

                if (
                    len(best_set) == set_size
                    and tuple(best_set + [score for score in best_report.values()]) not in selected_sets
                ):
                    selected_sets.add(tuple(best_set + [score for score in best_report.values()]))
                else:
                    available_oligos = adjusted_scores.sort_values(ascending=self.ascending).index.to_list()

            if len(selected_sets) == self.n_sets:
                break
        selected_sets = [list(selected_set) for selected_set in selected_sets]

        return selected_sets


class GraphBasedSelectionPolicy(OligoSelectionPolicy):
    def run_selection(self, oligos_scores: pd.Series, overlapping_matrix: csr_matrix, n_attempts: int):

        idx_to_oligo = {idx: oligo for idx, oligo in enumerate(oligos_scores.index)}
        G = nx.from_scipy_sparse_array(overlapping_matrix)
        G = nx.relabel_nodes(G, {i: idx_to_oligo[i] for i in G.nodes})

        # First check if there are no cliques with n oligos
        cliques = nx.algorithms.clique.find_cliques(G)
        n = self.opt_oligoset_size
        clique_init = []

        for clique in cliques:
            if len(clique) > self.min_oligoset_size:
                clique_init = clique
            if len(clique) >= n:
                break

        if not clique_init:
            # if no clique with min_oligoset_size was found we don't need to compute the sets
            return []

        n = min(n, len(clique_init))
        oligoset_init, oligoset_init_scores = self.set_scoring.apply(oligos_scores.loc[clique_init], n)

        # apply the heuristic
        clique_heuristic, oligos_scores = self._heuristic_selection(
            oligoset_init=oligoset_init,
            oligos_scores=oligos_scores,
            overlapping_matrix=overlapping_matrix,
            n_oligo=n,
            ascending=self.ascending,
        )
        overlapping_matrix_ids = list(oligos_scores.index)
        # overwrite initial oligoset
        oligoset_init, oligoset_init_scores = self.set_scoring.apply(oligos_scores.loc[clique_heuristic], n)
        # only keep oligos in overlap matrix that pass the heuristic
        overlapping_matrix_indices = [
            overlapping_matrix_ids.index(oligo_id) for oligo_id in oligos_scores.index
        ]
        overlapping_matrix = overlapping_matrix[overlapping_matrix_indices, :][:, overlapping_matrix_indices]
        overlapping_matrix_ids = [overlapping_matrix_ids[idx] for idx in overlapping_matrix_indices]

        # recompute graph and cliques from reduced matrix
        G = nx.from_scipy_sparse_array(overlapping_matrix)
        G = nx.relabel_nodes(G, {i: overlapping_matrix_ids[i] for i in range(len(overlapping_matrix_ids))})

        # need to recompute cliques to be able to reiterate through them from the start and find all sets
        cliques = nx.algorithms.clique.find_cliques(G)

        # Initialize oligoset results table
        oligosets = [list(oligoset_init) + list(oligoset_init_scores.values())]
        # Note: Search could be further optimised by iteratively throwing out oligos with worse scores then current best set
        for count, clique in enumerate(cliques):
            # Limit the number of combinations we iterate through
            if count > n_attempts:
                break
            if len(clique) >= n:
                # Get oligo_ids of clique, maybe create a function
                oligoset, oligoset_scores = self.set_scoring.apply(oligos_scores.loc[clique], n)
                oligosets.append(list(oligoset) + list(oligoset_scores.values()))

        return oligosets

    def _heuristic_selection(
        self,
        oligoset_init: list,
        oligos_scores: pd.Series,
        overlapping_matrix: csr_matrix,
        n_oligo: int,
        ascending: bool,
        n_trials: int = 100,
    ):

        # Sort the oligos by their score
        oligo_ids_sorted = oligos_scores.sort_values(ascending=ascending).index.to_list()

        # Sort overlap matrix by oligo scores, i.e. best performing oligo at first entry
        overlapping_matrix_ids = list(oligos_scores.index)
        overlapping_matrix_indices = [overlapping_matrix_ids.index(oligo_id) for oligo_id in oligo_ids_sorted]
        overlapping_matrix_sorted = overlapping_matrix[overlapping_matrix_indices, :][
            :, overlapping_matrix_indices
        ]

        # Initialize max_score with score from initial oligoset
        best_oligoset = oligoset_init
        best_oligoset_max_score = oligos_scores.loc[oligoset_init].max()

        for first_idx in range(min(len(oligo_ids_sorted), n_trials)):
            # Use the integer index because the matrix is converted to np array
            oligoset_idxs = np.array([first_idx])
            for _ in range(n_oligo - 1):
                # Find first oligo in sorted array that is not overlapping with any selected oligo
                no_overlap = np.all(overlapping_matrix_sorted[oligoset_idxs].toarray(), axis=0)
                # Ensure not to select already selected oligos
                # no_overlap[oligoset_idxs] = False
                # Add the first entry withou overlap to set, i.e. the next oligo with the best score
                if np.any(no_overlap):
                    oligoset_idxs = np.append(oligoset_idxs, np.where(no_overlap)[0][0])
                else:
                    break

            # If enough non overlapping oligos are found and score is lower than existing set, replace existing with new set
            if len(oligoset_idxs) == n_oligo:
                oligoset = [oligo_ids_sorted[idx] for idx in oligoset_idxs]
                oligoset_max_score = oligos_scores.loc[oligoset].max()
                if oligoset_max_score < best_oligoset_max_score:
                    best_oligoset_max_score = oligoset_max_score
                    best_oligoset = oligoset

        oligos_scores = oligos_scores[oligos_scores <= best_oligoset_max_score]

        return best_oligoset, oligos_scores