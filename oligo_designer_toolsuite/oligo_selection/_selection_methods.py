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
    """A class that defines the policy for selecting oligos based on a scoring system and pre-filtering criteria.

    The OligoSelectionPolicy class is responsible for selecting oligos based on the provided scoring function
    and optionally applying a pre-filtering step before the selection process.

    :param set_scoring: Scoring method for evaluating the quality of the oligo sets.
    :type set_scoring: SetScoringBase
    :param pre_filter: A flag indicating whether pre-filtering should be applied to the oligos before selection,
                which improves performance for larger sets (e.g., > 30) but can slow down small set selection (e.g., < 30).
    :type pre_filter: bool
    """

    def __init__(
        self,
        set_scoring: SetScoringBase,
        pre_filter: bool,
    ) -> None:
        """Constructor for the OligoSelectionPolicy class."""
        self.set_scoring = set_scoring
        self.pre_filter = pre_filter

    def apply(
        self,
        set_size_opt: int,
        set_size_min: int,
        n_sets: int,
        oligos_scores: pd.Series,
        non_overlap_matrix: csr_matrix,
        non_overlap_matrix_ids: list,
    ) -> pd.DataFrame:
        """
        Applies the selection policy to generate oligo sets based on the provided oligo scores and the non-overlap matrix.

        The method selects oligo sets from the given list of oligo ranked by their scores, applying pre-filtering if enabled to reduce the initial search space.
        Then, the optimal sets are selected using a custom selection algorithm, and the results are returned as a formatted DataFrame.

        :param set_size_opt: The optimal size of each oligo set.
        :type set_size_opt: int
        :param set_size_min: The minimum allowed size of each oligo set.
        :type set_size_min: int
        :param n_sets: The number of oligo sets to generate.
        :type n_sets: int
        :param oligos_scores: A pandas Series containing the scores for each oligo.
        :type oligos_scores: pd.Series
        :param non_overlap_matrix: A sparse matrix representing the inverse of the overlap between oligos.
        :type non_overlap_matrix: csr_matrix
        :param non_overlap_matrix_ids: A list of oligo IDs corresponding to the non-overlap matrix.
        :type non_overlap_matrix_ids: list
        :return: A DataFrame containing the selected oligo sets and their scores.
        :rtype: pd.DataFrame
        """
        self.opt_oligoset_size = set_size_opt
        self.min_oligoset_size = set_size_min
        self.n_sets = n_sets

        if self.pre_filter:
            oligos_scores, non_overlap_matrix, non_overlap_matrix_ids = self._pre_filter_oligos(
                oligos_scores, non_overlap_matrix, non_overlap_matrix_ids
            )

        if len(oligos_scores) > self.min_oligoset_size:
            oligosets = self._run_selection(oligos_scores, non_overlap_matrix, non_overlap_matrix_ids)
        else:
            oligosets = []

        # not clean but it works
        score_names = list(self.set_scoring.apply(pd.Series([], dtype=object), 0)[1].keys())

        return self._format_as_df(oligosets, score_names)

    def _pre_filter_oligos(
        self,
        oligos_scores: pd.Series,
        non_overlap_matrix: csr_matrix,
        non_overlap_matrix_ids: list,
    ) -> tuple[pd.Series, csr_matrix, list]:
        """
        Pre-filters oligos based on the minimal set removing all oligos from the initial set that
        form independent sets which are smaller than the minimal set size.

        :param oligos_scores: A pandas Series containing the scores for each oligo.
        :type oligos_scores: pd.Series
        :param non_overlap_matrix: A sparse matrix representing the inverse of the overlap between oligos.
        :type non_overlap_matrix: csr_matrix
        :param non_overlap_matrix_ids: A list of oligo IDs corresponding to the non-overlap matrix.
        :type non_overlap_matrix_ids: list
        :return: A tuple containing the filtered oligo scores and the updated overlap matrix and a list of oligo ids.
        :rtype: tuple[pd.Series, csr_matrix, list]
        """
        G = nx.from_scipy_sparse_array(non_overlap_matrix)
        max_clique = nx.approximation.max_clique(G)
        nodes = []
        while len(max_clique) > self.min_oligoset_size:
            nodes += list(max_clique)
            G.remove_nodes_from(max_clique)
            max_clique = nx.approximation.max_clique(G)

        # only keep oligos in overlap matrix that pass the pre-filter
        oligos_scores = oligos_scores.iloc[nodes]
        non_overlap_matrix_ids_tmp = [
            non_overlap_matrix_ids.index(oligo_id) for oligo_id in oligos_scores.index
        ]
        non_overlap_matrix = non_overlap_matrix[non_overlap_matrix_ids_tmp, :][:, non_overlap_matrix_ids_tmp]
        non_overlap_matrix_ids = [non_overlap_matrix_ids[idx] for idx in non_overlap_matrix_ids_tmp]
        return oligos_scores, non_overlap_matrix, non_overlap_matrix_ids

    @abstractmethod
    def _run_selection(
        self, oligos_scores: pd.Series, non_overlap_matrix: csr_matrix, non_overlap_matrix_ids: list
    ) -> list:
        """
        Selects the optimal oligo sets based on the provided oligo scores and overlap matrix.

        This is an abstract method that should be implemented by subclasses to define the specific selection algorithm.

        :param oligos_scores: A pandas Series containing the scores for each oligo.
        :type oligos_scores: pd.Series
        :param non_overlap_matrix: A sparse matrix representing the inverse of the overlap between oligos.
        :type non_overlap_matrix: csr_matrix
        :param non_overlap_matrix_ids: A list of oligo IDs corresponding to the non-overlap matrix.
        :type non_overlap_matrix_ids: list
        :return: A list of the selected oligo sets.
        :rtype: list
        """

    def _format_as_df(self, oligosets: list, score_names: list) -> pd.DataFrame:
        """
        Formats the oligo sets and scores into a pandas DataFrame.

        This method converts the list of oligo sets and their corresponding scores into a DataFrame, sorts the oligo sets by
        score, removes duplicates, and returns the top N oligo sets based on the selection criteria.

        :param oligosets: A list of selected oligo sets, each containing oligo IDs and their scores.
        :type oligosets: list
        :param score_names: A list of the names of the score columns to be used in the final DataFrame.
        :type score_names: list
        :return: A DataFrame containing the selected oligo sets and their scores, sorted by score.
        :rtype: pd.DataFrame
        """
        if len(oligosets) > 0:
            oligoset_size = len(oligosets[0]) - len(score_names)
            oligosets_columns = [f"oligo_{i}" for i in range(oligoset_size)] + score_names
            oligosets = pd.DataFrame(
                columns=oligosets_columns,
                data=oligosets,
            )

            # Sort oligosets by score
            oligosets.drop_duplicates(inplace=True, subset=oligosets.columns[:-1])
            oligosets.sort_values(
                list(oligosets.columns[oligoset_size:]),
                ascending=self.set_scoring.ascending,
                inplace=True,
            )
            oligosets = oligosets.head(self.n_sets)
            oligosets.reset_index(drop=True, inplace=True)
            oligosets.insert(0, "oligoset_id", oligosets.index)
            return oligosets
        else:
            return None


class GreedySelectionPolicy(OligoSelectionPolicy):
    """
    GreedySelectionPolicy selects oligo sets based on a greedy approach, where at each step, the best oligos are chosen
    according to a scoring function. A penalty can be applied to discourage the selection of overlapping oligos.

    The policy aims to select oligos that maximize the given scoring criteria while minimizing overlap between them.

    :param set_scoring: Scoring method for evaluating the quality of the oligo sets.
    :type set_scoring: SetScoringBase
    :param score_criteria: String representing the score metric to use for set selection (e.g., "set_score_worst", "set_score_sum", see SetScoringBase).
    :type score_criteria: str
    :param pre_filter: A flag indicating whether pre-filtering should be applied to the oligos before selection,
                which improves performance for larger sets (e.g., > 30) but can slow down small set selection (e.g., < 30).
    :type pre_filter: bool
    :param penalty: Penalty factor applied to selected oligos, defaults to 0.05 (the higher the value, the more distinct the sets).
    :type penalty: float, optional
    :param n_attempts: The number of attempts to make when generating oligo sets. Default is 1000.
    :type n_attempts: int, optional
    """

    def __init__(
        self,
        set_scoring: SetScoringBase,
        score_criteria: str,
        pre_filter: bool,
        penalty: float = 0.05,
        n_attempts: int = 1000,
    ) -> None:
        """Constructor for the GreedySelectionPolicy class."""
        super().__init__(set_scoring=set_scoring, pre_filter=pre_filter)
        self.score_criteria = score_criteria
        self.penalty = penalty
        self.n_attempts = n_attempts

    def _run_selection(
        self, oligos_scores: pd.Series, non_overlap_matrix: csr_matrix, non_overlap_matrix_ids: list
    ) -> list:
        """
        Selects a set of oligos based on a greedy algorithm that maximizes the scoring criteria while minimizing overlap
        between the selected oligo sets. The algorithm iteratively selects the best oligos for each set, adjusting the score
        for each oligo to encourage diversity and penalize overlap with previously selected oligos. It continues to select
        oligos until the desired number of sets is reached or the maximum number of attempts is made.

        :param oligos_scores: A pandas Series containing the scores for each oligo, where the index is the oligo identifier.
        :type oligos_scores: pd.Series
        :param non_overlap_matrix: A sparse matrix representing the inverse of the overlap between oligos.
        :type non_overlap_matrix: csr_matrix
        :param non_overlap_matrix_ids: A list of oligo IDs corresponding to the non-overlap matrix.
        :type non_overlap_matrix_ids: list
        :return: A list of lists, each containing a set of selected oligos.
        :rtype: list
        """
        oligo_to_idx = {oligo: idx for idx, oligo in enumerate(non_overlap_matrix_ids)}
        selected_sets = set()
        adjusted_scores = oligos_scores.copy()

        for set_size in range(self.opt_oligoset_size, self.min_oligoset_size - 1, -1):
            attempts = 0
            while len(selected_sets) < self.n_sets and attempts < self.n_attempts:
                attempts += 1
                best_set = []
                available_oligos = adjusted_scores.sort_values(
                    ascending=self.set_scoring.ascending
                ).index.to_list()
                for _ in range(set_size):
                    best_oligo = None
                    best_score = None

                    for oligo in available_oligos:
                        if not all(
                            non_overlap_matrix[oligo_to_idx[oligo], oligo_to_idx[best_set_oligo]]
                            for best_set_oligo in best_set
                        ):
                            continue

                        current_set = best_set + [oligo]
                        current_scores = adjusted_scores.loc[current_set]

                        _, scores = self.set_scoring.apply(current_scores, set_size)

                        score = scores[self.score_criteria]
                        better_than_best = (
                            best_score is None
                            or (self.set_scoring.ascending and score < best_score)
                            or (not self.set_scoring.ascending and score > best_score)
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
                            (1 + self.penalty) if self.set_scoring.ascending else (1 - self.penalty)
                        )

                if (
                    len(best_set) == set_size
                    and tuple(best_set + [score for score in best_report.values()]) not in selected_sets
                ):
                    selected_sets.add(tuple(best_set + [score for score in best_report.values()]))
                else:
                    available_oligos = adjusted_scores.sort_values(
                        ascending=self.set_scoring.ascending
                    ).index.to_list()

            if len(selected_sets) == self.n_sets:
                break
        selected_sets = [list(selected_set) for selected_set in selected_sets]

        return selected_sets


class GraphBasedSelectionPolicy(OligoSelectionPolicy):
    """
    Graph-Based Selection Policy for oligo selection. This policy utilizes a graph-based approach
    to select oligos, where the nodes represent the oligos and the edges represent their overlap wrt
    genomic coordinates. The selection process optimizes for a set of oligos that don't overlap and
    maximize the scoring criteria, while optionally using a heuristic approach to improve efficiency.

    :param set_scoring: Scoring method for evaluating the quality of the oligo sets.
    :type set_scoring: SetScoringBase
    :param pre_filter: A flag indicating whether pre-filtering should be applied to the oligos before selection,
                which improves performance for larger sets (e.g., > 30) but can slow down small set selection (e.g., < 30).
    :type pre_filter: bool
    :param n_attempts: The number of attempts to make when generating oligo sets. Default is 1000.
    :type n_attempts: int, optional
    :param heuristic: Whether to use a heuristic approach for faster selection. Defaults to True.
    :type heuristic: bool, optional
    :param heuristic_n_attempts: The number of attempts to find the optimal oligo set using the heuristic approach. Default is 1000.
    :type heuristic_n_attempts: int, optional
    :param clique_init_approximation: Whether to use an approximation approach for faster finding of initial oligo set. Defaults to False.
    :type clique_init_approximation: bool, optional
    """

    def __init__(
        self,
        set_scoring: SetScoringBase,
        pre_filter: bool,
        n_attempts: int = 1000,
        heuristic: bool = True,
        heuristic_n_attempts: int = 1000,
        clique_init_approximation=False,
    ) -> None:
        """Constructor for the GraphBasedSelectionPolicy class."""
        super().__init__(set_scoring=set_scoring, pre_filter=pre_filter)
        self.n_attempts = n_attempts
        self.heuristic = heuristic
        self.heuristic_n_attempts = heuristic_n_attempts
        self.clique_init_approximation = clique_init_approximation

    def _run_selection(
        self, oligos_scores: pd.Series, non_overlap_matrix: csr_matrix, non_overlap_matrix_ids: list
    ) -> list:
        """
        Runs the selection process for oligo sets based on the scores and overlap matrix. This method identifies cliques
        in the overlap graph, which represent non-overlapping sets of oligos wrt the their genomic coordinates, and then
        iteratively selects clique of oligos which contain at least n oligos. A heuristic can be applied to improve the
        efficiency of the graph-based selection.

        :param oligos_scores: A pandas Series containing the scores for each oligo, where the index is the oligo identifier.
        :type oligos_scores: pd.Series
        :param non_overlap_matrix: A sparse matrix representing the inverse of the overlap between oligos.
        :type non_overlap_matrix: csr_matrix
        :param non_overlap_matrix_ids: A list of oligo IDs corresponding to the non-overlap matrix.
        :type non_overlap_matrix_ids: list
        :return: A list of lists, each containing a set of selected oligos.
        :rtype: list
        """
        G = nx.from_scipy_sparse_array(non_overlap_matrix)
        G = nx.relabel_nodes(G, {i: non_overlap_matrix_ids[i] for i in range(len(non_overlap_matrix_ids))})

        # First check if there are no cliques with n oligos.
        # If we have large sets it's better to check the max clique size with an heuristic
        # because iterating through cliques until we find one which is greater than the
        # optimal set size, will take a long time.
        oligoset_size = self.opt_oligoset_size
        clique_init = []

        if self.clique_init_approximation:
            clique_max = nx.approximation.max_clique(G)
            if len(clique_max) > self.min_oligoset_size:
                clique_init = list(clique_max)
        else:
            cliques = nx.algorithms.clique.find_cliques(G)
            for clique in cliques:
                if len(clique) > self.min_oligoset_size:
                    clique_init = clique
                if len(clique) >= oligoset_size:
                    break

        if not clique_init:
            # if no clique with min_oligoset_size was found we don't need to compute the sets
            return []

        oligoset_size = min(oligoset_size, len(clique_init))
        oligoset_init, oligoset_init_scores = self.set_scoring.apply(
            oligos_scores.loc[clique_init], oligoset_size
        )

        if self.heuristic and oligoset_size == self.opt_oligoset_size:
            # apply the heuristic
            clique_heuristic, oligos_scores = self._heuristic_selection(
                oligoset_init=oligoset_init,
                oligos_scores=oligos_scores,
                non_overlap_matrix=non_overlap_matrix,
                non_overlap_matrix_ids=non_overlap_matrix_ids,
                oligoset_size=oligoset_size,
                heuristic_n_attempts=self.heuristic_n_attempts,
            )

            # overwrite initial oligoset
            oligoset_init, oligoset_init_scores = self.set_scoring.apply(
                oligos_scores.loc[clique_heuristic], oligoset_size
            )
            # only keep oligos in overlap matrix that pass the heuristic
            non_overlap_matrix_ids_tmp = [
                non_overlap_matrix_ids.index(oligo_id) for oligo_id in oligos_scores.index
            ]
            non_overlap_matrix = non_overlap_matrix[non_overlap_matrix_ids_tmp, :][
                :, non_overlap_matrix_ids_tmp
            ]
            non_overlap_matrix_ids = [non_overlap_matrix_ids[idx] for idx in non_overlap_matrix_ids_tmp]

            # recompute graph and cliques from reduced matrix
            G = nx.from_scipy_sparse_array(non_overlap_matrix)
            G = nx.relabel_nodes(
                G, {i: non_overlap_matrix_ids[i] for i in range(len(non_overlap_matrix_ids))}
            )

        # need to recompute cliques to be able to reiterate through them from the start and find all sets
        cliques = nx.algorithms.clique.find_cliques(G)

        # Initialize oligoset results table
        oligosets = [list(oligoset_init) + list(oligoset_init_scores.values())]
        # Note: Search could be further optimised by iteratively throwing out oligos with worse scores then current best set
        for count, clique in enumerate(cliques):
            # Limit the number of combinations we iterate through
            if count > self.n_attempts:
                break
            if len(clique) >= oligoset_size:
                # Get oligo_ids of clique, maybe create a function
                oligoset, oligoset_scores = self.set_scoring.apply(oligos_scores.loc[clique], oligoset_size)
                oligosets.append(list(oligoset) + list(oligoset_scores.values()))

        return oligosets

    def _heuristic_selection(
        self,
        oligoset_init: list,
        oligos_scores: pd.Series,
        non_overlap_matrix: csr_matrix,
        non_overlap_matrix_ids: list,
        oligoset_size: int,
        heuristic_n_attempts: int,
    ):
        """
        A heuristic approach to improve the selection of oligo sets by iteratively selecting non-overlapping oligos that maximize the
        score. This method selects the best set of oligos based on the score while ensuring that there is no overlap between the oligos.

        :param oligoset_init: The initial set of oligos to start the heuristic selection process.
        :type oligoset_init: list
        :param oligos_scores: A pandas Series containing the scores for each oligo, where the index is the oligo identifier.
        :type oligos_scores: pd.Series
        :param non_overlap_matrix: A sparse matrix representing the inverse of the overlap between oligos.
        :type non_overlap_matrix: csr_matrix
        :param non_overlap_matrix_ids: A list of oligo IDs corresponding to the non-overlap matrix.
        :type non_overlap_matrix_ids: list
        :param oligoset_size: The size of the oligo set to select.
        :type oligoset_size: int
        :param heuristic_n_attempts: The number of attempts to find the optimal oligo set using the heuristic approach.
        :type heuristic_n_attempts: int
        :return: The best set of oligos found by the heuristic approach and the updated oligos scores.
        :rtype: tuple(list, pd.Series)
        """
        # Sort the oligos by their score
        oligo_ids_sorted = oligos_scores.sort_values(ascending=self.set_scoring.ascending).index.to_list()

        # Sort overlap matrix by oligo scores, i.e. best performing oligo at first entry
        non_overlap_matrix_indices = [non_overlap_matrix_ids.index(oligo_id) for oligo_id in oligo_ids_sorted]
        non_overlap_matrix_sorted = non_overlap_matrix[non_overlap_matrix_indices, :][
            :, non_overlap_matrix_indices
        ]

        # Initialize max_score with score from initial oligoset
        best_oligoset = oligoset_init
        best_oligoset_max_score = oligos_scores.loc[oligoset_init].max()

        for first_idx in range(min(len(oligo_ids_sorted), heuristic_n_attempts)):
            # Use the integer index because the matrix is converted to np array
            oligoset_idxs = np.array([first_idx])
            for _ in range(oligoset_size - 1):
                # Find first oligo in sorted array that is not overlapping with any selected oligo
                no_overlap = np.all(non_overlap_matrix_sorted[oligoset_idxs].toarray(), axis=0)
                # Ensure not to select already selected oligos
                # no_overlap[oligoset_idxs] = False
                # Add the first entry withou overlap to set, i.e. the next oligo with the best score
                if np.any(no_overlap):
                    oligoset_idxs = np.append(oligoset_idxs, np.where(no_overlap)[0][0])
                else:
                    break

            # If enough non overlapping oligos are found and score is lower than existing set, replace existing with new set
            if len(oligoset_idxs) == oligoset_size:
                oligoset = [oligo_ids_sorted[idx] for idx in oligoset_idxs]
                oligoset_max_score = oligos_scores.loc[oligoset].max()
                if oligoset_max_score < best_oligoset_max_score:
                    best_oligoset_max_score = oligoset_max_score
                    best_oligoset = oligoset

        oligos_scores = oligos_scores[oligos_scores <= best_oligoset_max_score]

        return best_oligoset, oligos_scores
