############################################
# imports
############################################

import os
import re
import pandas as pd

from pathlib import Path
from joblib import Parallel, delayed
from Bio.Blast.Applications import NcbiblastnCommandline, NcbimakeblastdbCommandline

from . import SpecificityFilterBase, Blastn

import networkx as nx

############################################
# Oligo Blast Filter Classes
############################################

# TODO: remove inheritance


class CrossHybridizationFilter(SpecificityFilterBase):
    def __init__(
        self,
        dir_specificity: str,
        specificity_filter: SpecificityFilterBase,
        policy,
        n_jobs: int,
    ):
        """Constructor."""
        super().__init__(dir_specificity)

        self.dir_blast = os.path.join(self.dir_specificity, "blast")
        Path(self.dir_blast).mkdir(parents=True, exist_ok=True)

        self.dir_fasta = os.path.join(self.dir_specificity, "fasta")
        Path(self.dir_fasta).mkdir(parents=True, exist_ok=True)

        self.specificity_filter = specificity_filter
        self.n_jobs = n_jobs
        self.policy = policy

    def apply(self, oligo_database: dict, file_reference: str, n_jobs: int):
        oligos_per_region = self.count_oligos_per_region(oligo_database)
        cross_hybridization_graph = self.create_cross_hybridization_graph(
            oligo_database, file_reference, n_jobs=n_jobs
        )
        matching_oligos = self.policy(cross_hybridization_graph, oligos_per_region)

        # Can be parallelized if needed
        regions = list(oligo_database.keys())
        for region in regions:
            filtered_database_region = self._filter_matching_oligos(
                oligo_database[region], matching_oligos[region]
            )
            oligo_database[region] = filtered_database_region
        return oligo_database

    def create_cross_hybridization_graph(
        self, oligo_database: dict, reference_fasta: str, n_jobs: int
    ):
        database_name = self.specificity_filter.create_index(
            reference_fasta, n_jobs=n_jobs
        )
        matching_oligo_pairs = self.specificity_filter.get_all_matching_oligo_pairs(
            oligo_database, database_name, n_jobs
        )
        return nx.from_edgelist(matching_oligo_pairs)

    def count_oligos_per_region(self, oligo_database):
        return sorted(
            {
                region: len(oligo_database.database[region])
                for region in oligo_database.database.keys()
            }
        )


# Policies


def _oligo_to_region(oligo_name: str):
    return oligo_name.split("_")[0]


def remove_nodes_by_priority_policy(
    graph: nx.Graph, oligos_per_region: dict
) -> list(str):
    """removes oligos (graph nodes) based on a priorities dictionnary: nodes associated with regions with higher number oligos_per_region will be removed first. The process is repeated until there are no edges in the graph anymore.
    Args:
        graph (nx.Graph): graph representing oligos, an edge means that there is a match between two oligos
        oligos_per_region (dict): dictionnary that maps each region with the number of oligos
    Returns:
        list(str): list of oligos that are removed
    """

    removed_oligos = {region: [] for region in oligos_per_region.keys()}

    while graph.number_of_edges() > 0:
        edge = list(graph.edges)[0]
        region_0 = _oligo_to_region(edge[0])
        region_1 = _oligo_to_region(edge[1])
        if oligos_per_region[region_0] > oligos_per_region[region_1]:
            graph.remove_node(edge[0])
            removed_oligos[region_0].append(edge[0])
        else:
            graph.remove_node(edge[1])
            removed_oligos[region_1].append(edge[1])
    return removed_oligos


def remove_nodes_by_degree_policy(
    graph: nx.Graph, oligos_per_region: dict
) -> list(str):
    removed_oligos = {region: [] for region in oligos_per_region.keys()}

    while graph.number_of_edges() > 0:
        degrees = dict(graph.degree())

        if degrees:
            max_degree_node = max(degrees, key=degrees.get)
            graph.remove_node(max_degree_node)
            removed_oligos[_oligo_to_region(max_degree_node)].append(max_degree_node)
    return removed_oligos
