import os
from math import ceil
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd
from joblib import Parallel, delayed


class ProbesetGenerator:
    """This class is used to generate ranked, non-overlapping probe sets."""

    def __init__(
        self,
        n_probes_per_gene,
        min_n_probes_per_gene,
        probes_scoring,
        set_scoring,
        dir_probe_sets="probe_sets",  # write the output
        heurustic_selection=None,
    ) -> None:
        """Initialize the class.

        :param n_probes_per_gene: optimal number of probes that each set should contain
        :type n_probes_per_gene: int
        :param min_n_probes_per_gene: minimum number of probes that each set should contain
        :type min_n_probes_per_gene: int
        :param probes_scoring: class that scores hte probes and the sets of probes
        :type probes_scoring: ProbeScoring class
        :param dir_probe_sets: directory where the sets are written
        :type dir_probe_sets: str
        :param heurustic_selection: functions that preselects the probes
        :type dir_probe_sets: func
        """

        self.n_probes_per_gene = n_probes_per_gene
        self.min_n_probes_per_gene = min_n_probes_per_gene
        self.heurustic_selection = heurustic_selection
        self.probes_scoring = probes_scoring
        self.set_scoring = set_scoring
        self.dir_probe_sets = dir_probe_sets

    def get_probe_sets(self, DB, n_sets=50, n_jobs=None):
        """Generates in parallel the probesets and returns the DB class updated containing only the probes that belong to a set and with additional fields
        containing information about the sets, namely "set_id" and "set_score".

        :param DB: class containg the oligo sequences and their features
        :type DB: CustomDB class
        :param n_sets: maximal number of sets that will be generated, defaults to 50
        :type n_sets: int, optional
        :param n_jobs: nr of cores used, if None the value set in DB class is used, defaults to None
        :type n_jobs: int
        :return: Updated DB class
        :rtype: CustomDB class
        """

        # crete the folders where the files migh be written
        self.dir_probe_sets = os.path.join(DB.dir_output, self.dir_probe_sets)
        Path(self.dir_probe_sets).mkdir(parents=True, exist_ok=True)
        self.file_removed_genes = DB.file_removed_genes
        # set the number of cores
        if n_jobs is None:
            n_jobs = DB.n_jobs

        genes = list(DB.oligos_DB.keys())
        # generate batches
        genes_per_batch = ceil(len(genes) / n_jobs)
        genes_batches = [
            genes[genes_per_batch * i : min(len(genes) + 1, genes_per_batch * (i + 1))]
            for i in range(n_jobs)
        ]
        # get the probe set for this gene in parallel
        updated_oligos_DB = Parallel(n_jobs=n_jobs)(
            delayed(self._get_probe_set_for_batch)(batch, DB.oligos_DB, n_sets)
            for batch in genes_batches
        )
        # restore the oligo DB
        for batch, DB_batch in zip(genes_batches, updated_oligos_DB):
            for gene, probes in zip(batch, DB_batch):
                if probes is None:  # if some sets have been found
                    del DB.oligos_DB[gene]
                else:
                    DB.probesets[gene] = probes["probesets"]
                    del probes["probesets"]
                    DB.oligos_DB[gene] = probes

        """
        # get the probe set for this gene in parallel
        updated_oligos_DB = Parallel(n_jobs=n_jobs)( #there should be an explicit return
            delayed(self._get_probe_set_for_gene)(gene, DB.oligos_DB[gene], n_sets)
            for gene in genes
        )
        # restore the oligo DB
        for gene, probes in zip(genes, updated_oligos_DB):
            if probes is None:  # if some sets have been found
                del DB.oligos_DB[gene]
            else:
                DB.probesets[gene] = probes["probesets"]
                del probes["probesets"]
                DB.oligos_DB[gene] = probes
        """

        return DB

    def _get_probe_set_for_batch(self, batch, oligos_DB, n_sets):
        """Generate the probesets for the batch of genes. It returns a list of dictionaries containing the probes that belong to a set and with additional fields
        containing information about the sets, namely "set_id" and "set_score".

        :param batch: listo of the genes in teh batch
        :type batch: list
        :param oligos_DB: daatabse containing all the probes
        :type oligos_DB: class CustomDB
        :param n_sets: number of sets to generate
        :type n_sets: int
        :return: list of the updated dictionaries
        :rtype: list
        """
        updated_oligos_DB = []
        for gene in batch:
            updated_oligos_DB.append(
                self._get_probe_set_for_gene(gene, oligos_DB[gene], n_sets)
            )
        return updated_oligos_DB

    def _get_probe_set_for_gene(self, gene, probes, n_sets):
        """Generate the probesets for a gene. It returns the dictionary containing the probes that belong to a set and with additional fields
        containing information about the sets, namely "set_id" and "set_score".

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
            # gene is added to the lsit of genes with insufficient sets
            with open(self.file_removed_genes, "a") as handle:
                handle.write(f"{gene}\tOligo_selection\n")
            # gene is not required anymore
            return None
        else:
            probesets.to_csv(
                os.path.join(self.dir_probe_sets, f"probeset_{gene}.tsv"), sep="\t"
            )
            # update the dictionary adding the key sets for the probes in a set and  delleting the probes not included in any set
            updated_probes = {}
            for set_id, row in probesets.iterrows():
                for i in range(n):
                    if row[i] not in updated_probes:
                        updated_probes[row[i]] = probes[row[i]]
                        updated_probes[row[i]]["set_id"] = [set_id]
                        updated_probes[row[i]]["set_score"] = [list(row[n:])]
                    else:
                        updated_probes[row[i]]["set_id"].append(set_id)
                        updated_probes[row[i]]["set_score"].append(list(row[n:]))
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

        # probes_copy = copy.deepcopy(probes)  # probes is a mutable object (dictionary), hence is passed by reference
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
        n = self.n_probes_per_gene
        n_max = 0
        for clique in cliques:
            n_max = max(len(clique), n_max)
            if n_max >= n:
                break
        if n_max < n:
            if (
                n_max <= self.min_n_probes_per_gene
            ):  # in this case we don't need to compute the sets
                return n_max, None, None
            else:
                n = n_max
        # if we have an heuristic apply it
        heuristic_probeset = None
        if self.heurustic_selection is not None and n == self.n_probes_per_gene:
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

        return n, probesets, probes
