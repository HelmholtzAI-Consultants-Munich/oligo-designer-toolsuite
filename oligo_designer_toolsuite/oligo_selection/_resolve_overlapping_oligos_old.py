import os
from abc import ABC, abstractmethod
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd
from joblib import Parallel, delayed


# TODO
# - how are the data of the previous step passed? file or dictionary?
# - should we keep also the data frames containig the probesets? In case how can we store them unsing joblib?
class ProbesetGenerator:
    """This class is used to generate ranked, non-overlapping probe sets."""

    def __init__(
        self,
        n_probes_per_gene,
        min_n_probes_per_gene,
        probes_scoring,
        n_jobs=2,
        dir_output="output",
        dir_probes="probes",
        dir_probe_sets="probe_sets",
        heurustic_selection=None,
    ) -> None:
        """Initialize the class. The Scoring of the sets and the heuristic function is"""

        self.n_probes_per_gene = n_probes_per_gene
        self.min_n_probes_per_gene = min_n_probes_per_gene
        self.n_jobs = n_jobs
        self.dir_probes = os.path.join(dir_output, dir_probes)
        self.dir_probe_sets = os.path.join(dir_output, dir_probe_sets)
        Path(self.dir_probe_sets).mkdir(parents=True, exist_ok=True)
        self.file_removed_genes = os.path.join(
            dir_output, "genes_with_insufficient_probes.txt"
        )
        self.heurustic_selection = heurustic_selection
        self.probes_scoring = probes_scoring
        self.time = 0

    def get_probe_sets(self, n_sets=50):
        """Generates the sets for all the genes in parallel, for each gene we create a process"""
        # 1. Get the list of genes adn their probes
        # can we assume the we have a single folder for teh probes and each gene as an unique file?
        probes_files = [
            f for f in os.listdir(self.dir_probes) if f.startswith("probes_")
        ]
        genes_list = []
        probes_list = []
        for probe_file in probes_files:
            genes_list.append(probe_file.split("probes_", -1)[1].split(".")[0])
            probes_list.append(
                pd.read_csv(
                    os.path.join(self.dir_probes, probe_file), sep="\t", index_col=0
                )
            )
        # get the probe set for this gene in parallel
        Parallel(n_jobs=self.n_jobs)(
            delayed(self._get_probe_set_for_gene)(gene, probes, n_sets)
            for gene, probes in zip(genes_list, probes_list)
        )

    def _get_probe_set_for_gene(self, gene, probes, n_sets):
        """Create the sets for a specific gene, firstly the overlapping matrix is created and then the
        non overlapping sets are ranked and saved"""
        # create the overlapping matrix
        overlapping_matrix = self._get_overlapping_matrix(probes)
        # create the set
        n, probesets = self._get_non_overlapping_sets(
            probes, overlapping_matrix, n_sets
        )
        # delete the useless variable to free some memory(overlapping matrix)
        del overlapping_matrix  # free some memory
        # write the set
        if probesets is None:  # value passed as a parameter
            # gene is added to the lsit of genes with insufficient sets
            # TODO: write the genes into a file (should we include also the length)
            print(f"{gene} has insufficient sets")
            with open(self.file_removed_genes, "a") as handle:
                handle.write(f"{gene}\t{n}\n")
        else:
            probesets.to_csv(
                os.path.join(self.dir_probe_sets, f"probeset_{gene}.tsv"), sep="\t"
            )

    def _get_overlapping_matrix(self, probes):
        """Generate overlap matrix for the given gene"""

        def _get_overlap(seq1_intervals, seq2_intervals):
            for a in seq1_intervals:
                for b in seq2_intervals:
                    if min(a[1], b[1]) - max(a[0], b[0]) > -1:
                        return True
            return False

        probes_copy = probes.astype(
            "string"
        )  # probes is a mutable object (DataFrame), hence is passed by reference (can we avoid to copy the whole df?)
        probes_copy["starts"] = probes_copy["start"].str.split(";")
        probes_copy["ends"] = probes_copy["end"].str.split(";")
        probes_copy["intervals"] = [
            [[int(start[i]), int(end[i])] for i in range(len(start))]
            for start, end in zip(probes_copy["starts"], probes_copy["ends"])
        ]  # save a list of couples of [start,end] of the duplicates of that probe

        intervals = probes_copy.intervals.values

        matrix = np.zeros(
            (len(intervals), len(intervals)), dtype=int
        )  # on the diagonal we have overlaps
        for i in range(len(intervals)):
            for j in range(i + 1, len(intervals)):
                if _get_overlap(intervals[i], intervals[j]):
                    matrix[i, j] = 1
        matrix = matrix + np.transpose(matrix) + np.eye(len(intervals))
        matrix = np.ones((len(intervals), len(intervals)), dtype=int) - matrix

        matrix = pd.DataFrame(
            matrix, columns=probes.index, index=probes.index
        )  # probeid is already the index of the matrix
        return matrix

    def _get_non_overlapping_sets(self, probes, overlapping_matrix, n_sets):
        n = self.n_probes_per_gene
        # Score probes
        # review the application of the scoring class
        probes = self.probes_scoring.apply(probes)  # add a column score to the probes

        # Represent overlap matrix as graph
        GC = nx.convert_matrix.from_numpy_matrix(overlapping_matrix.values)
        # GC = nx.algorithms.operators.unary.complement(G)

        # Initialize variable
        heuristic_probeset = None

        # First check if there are no cliques with n probes
        cliques = nx.algorithms.clique.find_cliques(GC)
        n_max = 0
        for clique in cliques:
            n_max = max(len(clique), n_max)
            if n_max >= n:
                break
        if n_max < n:
            if (
                n_max < self.min_n_probes_per_gene
            ):  # in this case we don't need to compute the sets
                return n_max, None
            else:
                n = n_max
        # if we have an heuristic apply it
        if (
            self.heurustic_selection is not None and n == self.n_probes_per_gene
        ):  # apply the heuristic
            # heuristic should return already the filtered df
            probes, heuristic_probeset = self.heurustic_selection(
                probes, overlapping_matrix, n
            )
            heuristic_probeset = self.probes_scoring.get_probeset(
                heuristic_probeset, n
            )  # make it a list as for all the other cliques for future use
            # recompute the cliques
            GC = nx.convert_matrix.from_numpy_matrix(
                overlapping_matrix.loc[probes.index, probes.index].values
            )
            # GC = nx.algorithms.operators.unary.complement(G)
        cliques = nx.algorithms.clique.find_cliques(GC)

        probesets = []
        tot = 0
        # Note: Search could be further optimised by iteratively throwing out probes with worse scores then current best set
        for count, clique in enumerate(cliques):
            # Limit the number of combinations we iterate through
            if count > 100000:  # set a meaningful value
                break
            if len(clique) >= n:
                # Get probe_ids of clique
                tot += 1
                clique_probes = probes.iloc[clique]
                probeset = self.probes_scoring.get_probeset(clique_probes, n)
                probesets.append(probeset)
        # put the sets in a dataframe
        probesets = pd.DataFrame(
            columns=[f"probe_{i}" for i in range(n)] + ["set_score"], data=probesets
        )

        # add the heurustuc best set, if is in the best n-sets then it will be kept
        if heuristic_probeset:
            probesets.loc[len(probesets)] = heuristic_probeset
        # Sort probesets by score
        probesets = probesets.sort_values("set_score", ascending=True)
        probesets = probesets.head(n_sets)
        probesets = probesets.reset_index(drop=True)

        return n, probesets


# Is this scheme generalizable (each probe has a score and the score set is obtaind form the score of each sequence)?
class ProbeScoring(ABC):
    def __init__(self):
        pass

    def apply(self, probes):
        probes["probe_score"] = probes.apply(self.scoring_function, axis=1)
        return probes

    @abstractmethod
    def scoring_function(self, probe):
        # given a clique, creates the set from teh clique and returns the score/scores of the set
        pass

    @abstractmethod
    def get_probeset(clique_probes, n):
        # scores a set of probes, need to to have the score defined, otherwise it is computes
        pass


class PadlockScoring(ProbeScoring):
    def __init__(
        self, Tm_min, Tm_opt, Tm_max, GC_min, GC_opt, GC_max, Tm_weight=1, GC_weight=1
    ):
        # pass parameters as a list?
        super().__init__()
        self.Tm_min = Tm_min
        self.Tm_opt = Tm_opt
        self.Tm_max = Tm_max
        self.GC_min = GC_min
        self.GC_opt = GC_opt
        self.GC_max = GC_max
        self.Tm_weight = Tm_weight
        self.GC_weight = GC_weight
        self.generate_scoring_functions()

    def scoring_function(self, probe):
        # distance from the optimal melting temperature weightend by the how far is the optimum from the min/ max
        # teh scoring is the lower teh better
        Tm_dif = (
            probe["melting_temperature"] - self.Tm_opt
        )  # check the names of the columns
        GC_dif = probe["GC_content"] - self.GC_opt
        return self.Tm_weight * self.Tm_error(Tm_dif) + self.GC_weight * self.GC_error(
            GC_dif
        )

    def generate_scoring_functions(self):
        # in a symmetric case can become more efficint to evaluate
        # define the error function for the melting temperature
        Tm_dif_max = self.Tm_max - self.Tm_opt
        Tm_dif_min = self.Tm_opt - self.Tm_min
        if Tm_dif_max == Tm_dif_min:
            self.Tm_error = lambda Tm_dif: abs(Tm_dif) / Tm_dif_max
        else:
            self.Tm_error = lambda Tm_dif: abs(Tm_dif) / Tm_dif_max * (
                Tm_dif > 0
            ) + abs(Tm_dif) / Tm_dif_min * (Tm_dif < 0)
        # define the error function for the GC content
        GC_dif_max = self.GC_max - self.GC_opt
        GC_dif_min = self.GC_opt - self.GC_min
        if GC_dif_max == GC_dif_min:
            self.GC_error = lambda GC_dif: abs(GC_dif) / GC_dif_max
        else:
            self.GC_error = lambda GC_dif: abs(GC_dif) / GC_dif_max * (
                GC_dif > 0
            ) + abs(GC_dif) / GC_dif_min * (GC_dif < 0)

    def get_probeset(self, clique_probes, n):
        # retunrs a list with n+1 elements, first n are the pid of the probes, the last is the score of the set
        # the socre of teh set cal be a lisof scores to avoid ties in teh sorting. However the fosrt element should be the most relevant and are checked one by one
        # lists are comapred in lexographcal order, first compare the fisrt element, if they are equal, compare the second element and so on
        # check that we have the socres defined
        best_n_probes = clique_probes.sort_values("probe_score", ascending=True).head(n)

        # Calculate performance of probeset
        probeset_error_max = best_n_probes["probe_score"].max()
        probeset_error_sum = best_n_probes["probe_score"].sum()
        probeset = best_n_probes.index.tolist()
        probeset.append([probeset_error_max, probeset_error_sum])
        return probeset


def padlock_heuristic_selection(probes, overlapping_matrix, n_probe, n_trials=10000):
    # compute a heuristic selection and return the best set and the filtered probeset
    # returns also the best probeset found with the heuristic (df containing those sequences)

    # assert that we have a score column?
    probes_sorted = probes.sort_values(
        "probe_score"
    )  # used only to compute the best set
    probe_ids_sorted = probes_sorted.index.tolist()

    inv_mat_sorted = overlapping_matrix.loc[
        probe_ids_sorted, probe_ids_sorted
    ].values  # already done before  in teh code, should compute the complementary here again?
    # inv_mat_sorted = np.abs(inv_mat_sorted - 1)

    max_score = (
        probes_sorted.iloc[-1]["probe_score"] * 1.1
    )  # already sorted df, the max is the last entry
    best_idx_set = []
    for first_idx in range(
        min(len(probe_ids_sorted), n_trials)
    ):  # use the integer index because the matric is a np array
        set_idxs = np.array([first_idx])
        for _ in range(n_probe - 1):
            # find first probe in sorted array that is not overlapping with any selected probe
            no_overlap = np.all(inv_mat_sorted[set_idxs], axis=0)
            if np.any(no_overlap):
                set_idxs = np.append(set_idxs, np.where(no_overlap)[0][0])
            else:
                break
        if len(set_idxs) == n_probe:
            score = np.max(probes_sorted["probe_score"].values[set_idxs])
            if score < max_score:
                max_score = score
                best_idx_set = set_idxs
    best_set = probes_sorted.iloc[best_idx_set]

    probes_below_err = probes["probe_score"] <= max_score
    probes = probes.loc[probes_below_err]
    return probes, best_set
