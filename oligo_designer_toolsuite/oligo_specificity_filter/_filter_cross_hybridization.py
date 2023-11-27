############################################
# imports
############################################

import os
from pathlib import Path

import networkx as nx
from joblib import Parallel, delayed

from . import AlignmentSpecificityFilter, SpecificityFilterBase

############################################
# Oligo Blast Filter Classes
############################################


class CrossHybridizationFilter(SpecificityFilterBase):
    """
    This class implements a cross-hybridization filter for oligos. It constructs a graph representation of oligo interactions using an alignment method, then applies a policy to remove oligos which cross-hybridize.
    The class provides methods to create a cross-hybridization graph, and apply the chosen policy to filter out oligos.
    The policy is a function that determines which oligos to remove in case of cross-hybridization.

    :param specificity_filter: An instance of AlignmentSpecificityFilter class that provides matching oligo pairs based on an alignment criteria.
    :type specificity_filter: AlignmentSpecificityFilter
    :param policy: A function that defines the strategy for oligo removal. It takes the graph and a dictionary where keys are the databse regions as inputs and returns a dictionary of removed oligos.
    :type policy: function
    :param dir_cross_hybridization: Directory where files produced during the cross hybridization check are stored.
    :type dir_cross_hybridization: str
    :param n_jobs: Number of simultaneous parallel computations.
    :type n_jobs: int
    """

    def __init__(
        self,
        specificity_filter: AlignmentSpecificityFilter,
        policy,
        dir_cross_hybridization,
        n_jobs: int,
    ):
        self.specificity_filter = specificity_filter
        self.n_jobs = n_jobs
        self.policy = policy
        self.dir_cross_hybridization = dir_cross_hybridization

    def apply(self, database: dict, file_reference: str, n_jobs: int):
        oligos_per_set = self._count_oligos_per_set(database)
        regions = list(database.keys())

        # Process in parallel per region
        matching_oligos = Parallel(n_jobs=n_jobs)(
            delayed(self._apply_single_region)(database[region], region, oligos_per_set)
            for region in regions
        )

        for region in regions:
            filtered_database_region = self._filter_matching_oligos(
                database[region],
                matching_oligos[region],
            )
            database[region] = filtered_database_region
        return database

    def _create_cross_hybridization_graph(
        self, oligo_database: dict, reference_fasta: str, region: str
    ):
        database_name = self.specificity_filter.create_index(reference_fasta, n_jobs=1)
        matching_oligo_pairs = self.specificity_filter.get_matching_oligo_pairs(
            oligo_database, database_name, region
        )
        return nx.from_edgelist(matching_oligo_pairs)

    def _count_oligos_per_set(self, oligo_database):
        return {
            region: len(oligo_database.database[region])
            for region in oligo_database.database.keys()
        }

    def _apply_single_region(
        self,
        database,
        region,
        oligos_per_set,
    ):
        file_reference = database.write_fasta_from_database(
            self.dir_cross_hybridization,
            f"oligo_database_fasta_{region}",
            region=region,
        )
        cross_hybridization_graph = self._create_cross_hybridization_graph(
            database, file_reference, region
        )
        matching_oligos = self.policy(cross_hybridization_graph, oligos_per_set)
        return matching_oligos


# Policies


def remove_nodes_by_priority_policy(graph: nx.Graph, oligos_per_set: dict) -> list(str):
    """
    Removes oligos (graph nodes) based on a priorities dictionnary. The process is repeated until there are no edges in the graph anymore.
    :param graph: Graph where nodes represent oligos and edges indicate cross-hybridization between pairs of oligos.
    :type graph: nx.Graph
    :param oligos_per_set: Dictionary mapping each set to a priority. Nodes associated with sets with higher priority will be removed first
    :type oligos_per_set: dict
    :returns: Dictionary where each key is a set and the value is a list of oligos removed from that set.
    :rtype: dict[str, list[str]]
    """

    removed_oligos = {oligo_set: [] for oligo_set in oligos_per_set.keys()}

    while graph.number_of_edges() > 0:
        edge = list(graph.edges)[0]
        set_0 = _oligo_to_set(edge[0])
        set_1 = _oligo_to_set(edge[1])
        if oligos_per_set[set_0] > oligos_per_set[set_1]:
            graph.remove_node(edge[0])
            removed_oligos[set_0].append(edge[0])
        else:
            graph.remove_node(edge[1])
            removed_oligos[set_1].append(edge[1])
    return removed_oligos


def remove_nodes_by_degree_policy(graph: nx.Graph, oligos_per_set: dict) -> list(str):
    """
    Iteratively removes nodes with the highest degree from a graph until there are no edges in the graph anymore.

    :param graph: Graph where nodes represent oligos and edges indicate cross-hybridization between pairs of oligos.
    :type graph: nx.Graph
    :param oligos_per_set: Dictionary mapping each set with an arbitrary value (we use this structure to romain consistant with other functions).
    :type oligos_per_set: dict
    :returns: Dictionary where each key is a set and the value is a list of oligos removed from that set.
    :rtype: dict[str, list[str]]
    """
    removed_oligos = {oligo_set: [] for oligo_set in oligos_per_set.keys()}

    while graph.number_of_edges() > 0:
        degrees = dict(graph.degree())

        if degrees:
            max_degree_node = max(degrees, key=degrees.get)
            graph.remove_node(max_degree_node)
            removed_oligos[_oligo_to_set(max_degree_node)].append(max_degree_node)
    return removed_oligos
