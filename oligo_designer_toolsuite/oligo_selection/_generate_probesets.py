from typing import Callable

import networkx as nx
import numpy as np
import pandas as pd
from joblib import Parallel, delayed

from ..IO import CustomOligoDB
from ..oligo_efficiency import ProbeScoringBase, SetScoringBase


class ProbesetGenerator:
    """This class is used to generate ranked, non-overlapping probe sets. Each probe is scored according to the given ``probes_scoring`` class
    and the sets according to the given ``set_scoring`` calss.

    :param probeset_size: optimal number of probes that each set should contain
    :type probeset_size: int
    :param min_probeset_size: minimum number of probes that each set should contain
    :type min_probeset_size: int
    :param probes_scoring: class that scores hte probes and the sets of probes
    :type probes_scoring: ProbeScoring class
    :param dir_probe_sets: directory where the sets are written
    :type dir_probe_sets: str
    :param heurustic_selection: functions that preselects the probes making the seach of the best probesets less demanding
    :type dir_probe_sets: Callable
    :param write_genes_with_insufficient_probes: if True genes with insufficient probes are written in a file, defaults to True
    :type write_genes_with_insufficient_probes: bool, optional"""

    def __init__(
        self,
        probeset_size: int,
        min_probeset_size: int,
        probes_scoring: ProbeScoringBase,
        set_scoring: SetScoringBase,
        heurustic_selection: Callable = None,
        write_genes_with_insufficient_probes: bool = True,
    ) -> None:
        """Initialize the class."""

        self.probeset_size = probeset_size
        self.min_probeset_size = min_probeset_size
        self.heurustic_selection = heurustic_selection
        self.probes_scoring = probes_scoring
        self.set_scoring = set_scoring
        self.write_genes_with_insufficient_probes = write_genes_with_insufficient_probes

    def apply(self, database: CustomOligoDB, n_sets: int = 50, n_jobs: int = None):
        """Generates in parallel the probesets and selects the best ``n_sets`` according to the
        The database class is updated, in particular form the ``oligos_DB`` are filtered out all the probes that don't belong to any probeset and in the class attruibute ``probesets`` are stored
        the computed probesets. The latter is a dictionary having as keys the genes names and as values a pandas.DataFrame containinig the probesets. The strucutre of the pandas.DataFrame is the following:

        +-------------+----------+----------+----------+-------+----------+-------------+-------------+-------+
        | probeset_id | probe_0  | probe_1  | probe_2  |  ...  | probe_n  | set_score_1 | set_score_2 |  ...  |
        +-------------+----------+----------+----------+-------+----------+-------------+-------------+-------+
        | 0           | AGRN_184 | AGRN_133 | AGRN_832 |  ...  | AGRN_706 | 0.3445      | 1.2332      |  ...  |
        +-------------+----------+----------+-----+----+-------+----------+-------------+-------------+-------+

        :param database: class containg the oligo sequences and their features
        :type database: CustomDB class
        :param n_sets: maximal number of sets that will be generated, defaults to 50
        :type n_sets: int, optional
        :param n_jobs: nr of cores used, if None the value set in database class is used, defaults to None
        :type n_jobs: int
        :return: Updated database class
        :rtype: CustomDB class
        """

        # set the number of cores
        if n_jobs is None:
            n_jobs = database.n_jobs

        genes = list(database.oligos_DB.keys())
        # get the probe set for this gene in parallel
        updated_oligos_DB = Parallel(
            n_jobs=n_jobs
        )(  # there should be an explicit return
            delayed(self._get_probe_set_for_gene)(
                gene, database.oligos_DB[gene], n_sets
            )
            for gene in genes
        )
        # restore the oligo database
        for gene, probes in zip(genes, updated_oligos_DB):
            if probes is None:  # if some sets have been found
                database.oligos_DB[gene] = {}  # probeset is not generated
            else:
                database.probesets[gene] = probes["probesets"]
                del probes["probesets"]
                database.oligos_DB[gene] = probes
        database.remove_genes_with_insufficient_probes(
            pipeline_step="probeset generation",
            write=self.write_genes_with_insufficient_probes,
        )

        return database

    def _get_probe_set_for_gene(self, gene: str, probes: dict, n_sets: int):
        """Generate the probesets for a gene.

        :param gene: gene for whihc the probesets are computed
        :type gene: str
        :param probes: dictionary with the information about the probes
        :type probes: dict
        :param n_sets: number of sets to generate
        :type n_sets: int
        :return: updated_probes
        :rtype: dict
        """
        # create the overlapping matrix
        overlapping_matrix = self._get_overlapping_matrix(probes)
        # create the set
        n, probesets, probes = self._get_non_overlapping_sets(
            probes, overlapping_matrix, n_sets
        )
        # delete the useless variable to free some memory(overlapping matrix)
        del overlapping_matrix  # free some memory
        # write the set
        if probesets is None:  # value passed as a parameter
            return None  # no more probes are left for this gene
        else:
            # update the dictionary adding the key sets for the probes in a set and  delleting the probes not included in any set
            updated_probes = {}
            for set_id, row in probesets.iterrows():
                for i in range(1, n + 1):
                    if row[i] not in updated_probes:
                        updated_probes[row[i]] = probes[row[i]]
            # add also the probesets
            updated_probes["probesets"] = probesets
            return updated_probes

    def _get_overlapping_matrix(self, probes):
        """Creates a matrix that encodes the overlapping of different probes. the matrix has dimensions n_probes * n_probes where each
        row and column belong to a probe. Each entry contains 1 if the the correspondent probes don't overlap and 0 if they overlap, this
        is done because in the next this matrix will be used as an adjacency matrix and the sets of non overlapping probes are cliques of this graph.

        :pram probes: dictionary containing all the probes
        :type probes: dict
        :return: overlapping matrix
        :rtype: pandas.DataFrame
        """

        def _get_overlap(seq1_intervals, seq2_intervals):
            for a in seq1_intervals:
                for b in seq2_intervals:
                    if min(a[1], b[1]) - max(a[0], b[0]) > -1:
                        return True
            return False

        intervals = []
        probes_indices = []
        for probe_id in probes.keys():
            probes_indices.append(probe_id)  # keep track of the indices
            interval = []
            for start, end in zip(probes[probe_id]["start"], probes[probe_id]["end"]):
                interval.append(
                    [start, end]
                )  # save a list of couples of [start,end] of the duplicates of that probe
            intervals.append(interval)

        overlapping_matrix = np.zeros(
            (len(intervals), len(intervals)), dtype=int
        )  # on the diagonal we have overlaps
        for i in range(len(intervals)):
            for j in range(i + 1, len(intervals)):
                if _get_overlap(intervals[i], intervals[j]):
                    overlapping_matrix[i, j] = 1
        overlapping_matrix = (
            overlapping_matrix
            + np.transpose(overlapping_matrix)
            + np.eye(len(intervals))
        )
        overlapping_matrix = (
            np.ones((len(intervals), len(intervals)), dtype=int) - overlapping_matrix
        )
        overlapping_matrix = pd.DataFrame(
            data=overlapping_matrix,
            columns=probes_indices,
            index=probes_indices,
            dtype=int,
        )

        return overlapping_matrix

    def _get_non_overlapping_sets(self, probes, overlapping_matrix, n_sets):
        """Generates the non overlapping sets and return the best n_sets. Firstly the probes are scored and then, if it is available,
        an heuristic method is used to reduce the number of probes to the more promising ones. Then all the possible combination of non overlapping
        sets are considered and the best n-sets are returned.

        :param probes: dictionary containing all the probes
        :type probes: dict
        :param overlapping_matrix: matrix containig information about which probes overlap
        :type overlapping_matrix: pandas.DataFrame
        :param n_sets: number of sets to generate
        :type n_sets: int
        :return: actual size of the sets, DataFrame containnig the computed probesets, dictionary of probes
        :rtype: int, pandas.DataFrame, dict
        """
        probe_indices = np.array(overlapping_matrix.index)
        # Score probes and create a pd series with the same order of probes as in overlapping matrix
        probes, probes_scores = self.probes_scoring.apply(
            probes, probe_indices
        )  # add a entry score to the probes
        # Represent overlap matrix as graph
        G = nx.convert_matrix.from_numpy_array(overlapping_matrix.values)
        # First check if there are no cliques with n probes
        cliques = nx.algorithms.clique.find_cliques(G)
        n = self.probeset_size
        n_max = 0
        for clique in cliques:
            n_max = max(len(clique), n_max)
            if n_max >= n:
                break
        if n_max < n:
            if (
                n_max <= self.min_probeset_size
            ):  # in this case we don't need to compute the sets
                return n_max, None, None
            else:
                n = n_max
        # if we have an heuristic apply it
        heuristic_probeset = None
        if self.heurustic_selection is not None and n == self.probeset_size:
            # apply the heuristic
            probes, probes_scores, heuristic_set = self.heurustic_selection(
                probes, probes_scores, overlapping_matrix, n
            )
            heuristic_probeset = self.set_scoring.apply(
                heuristic_set, n
            )  # make it a list as for all the other cliques for future use
            # recompute the cliques
            G = nx.convert_matrix.from_numpy_array(
                overlapping_matrix.loc[probes_scores.index, probes_scores.index].values
            )
        # recompute the cliques
        cliques = nx.algorithms.clique.find_cliques(G)

        # Search the best set
        probesets = []
        # Note: Search could be further optimised by iteratively throwing out probes with worse scores then current best set
        for count, clique in enumerate(cliques):
            # Limit the number of combinations we iterate through
            if count > 100000:  # set a meaningful value
                break
            if len(clique) >= n:
                # Get probe_ids of clique, maybe create a function
                clique_probes = probes_scores.iloc[clique]
                probeset = self.set_scoring.apply(clique_probes, n)
                probesets.append(probeset)
        # put the sets in a dataframe
        probesets = pd.DataFrame(
            columns=[f"probe_{i}" for i in range(n)]
            + [f"set_score_{i}" for i in range(len(probesets[0]) - n)],
            data=probesets,
        )
        # add the heurustuc best set, if is in the best n-sets then it will be kept
        if heuristic_probeset:
            probesets.loc[len(probesets)] = heuristic_probeset
        # Sort probesets by score
        probesets.drop_duplicates(inplace=True, subset=probesets.columns[:-1])
        probesets.sort_values(list(probesets.columns[n:]), ascending=True, inplace=True)
        probesets = probesets.head(n_sets)
        probesets.reset_index(drop=True, inplace=True)
        probesets.insert(0, "probeset_id", probesets.index)

        return n, probesets, probes
