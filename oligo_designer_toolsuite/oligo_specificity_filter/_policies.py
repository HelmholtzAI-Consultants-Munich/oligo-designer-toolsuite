############################################
# imports
############################################

from abc import ABC, abstractmethod

import pandas as pd
import networkx as nx

from oligo_designer_toolsuite._constants import SEPARATOR_OLIGO_ID
from oligo_designer_toolsuite.database import OligoDatabase

############################################
# Hybridization Policies
############################################


class FilterPolicyBase(ABC):
    """An abstract base class for defining the policy for managing hybridization within an oligonucleotide database.
    Implementations of this class should provide specific strategies to minimize hybridization risk in
    the design of oligonucleotide sequences.
    """

    def __init__(self):
        """Constructor for the FilterPolicyBase class."""

    @abstractmethod
    def apply(self, graph: nx.Graph, oligo_database: OligoDatabase):
        """Applies the hybridization policy to an oligonucleotide database, removing oligos to reduce
        hybridization risks.

        :param graph: A networkx graph where nodes represent oligonucleotides and edges potential hybridizations.
        :type graph: nx.Graph
        :param oligo_database: The database of oligonucleotides to be evaluated.
        :type oligo_database: OligoDatabase
        """

    def _oligo_to_region(self, oligo: str) -> str:
        """Extracts the region identifier from an oligonucleotide ID.

        :param oligo: The oligonucleotide ID.
        :type oligo: str
        :return: The region identifier.
        :rtype: str
        """
        return oligo.split(SEPARATOR_OLIGO_ID)[0]

    def _get_number_oligos_per_region(self, oligo_database: OligoDatabase):
        """Calculates the number of oligonucleotides per region within the oligo database.

        :param oligo_database: The database of oligonucleotides.
        :type oligo_database: OligoDatabase
        :return: A dictionary mapping each region to its count of oligonucleotides.
        :rtype: dict
        """
        return {region: len(oligo_database.database[region]) for region in oligo_database.database.keys()}


class RemoveAllPolicy(FilterPolicyBase):
    """A policy that removes all oligonucleotides that potentially hybridize from the database,
    targeting both query and reference oligos identified in oligo pair hits.

    This policy aggressively ensures the elimination of any oligonucleotide involved in potential hybridization,
    removing them from their respective regions in the database.
    """

    def __init__(self):
        """Constructor for the RemoveAllPolicy class."""

    def apply(self, oligo_pair_hits: pd.DataFrame, oligo_database: OligoDatabase) -> dict[str, list[str]]:
        """Applies the removal policy to the oligonucleotide database based on the pair hits data frame.

        :param oligo_pair_hits: Data frame containing pairs of oligonucleotides that potentially hybridize.
        :type oligo_pair_hits: pd.DataFrame
        :param oligo_database: The database of oligonucleotides.
        :type oligo_database: OligoDatabase
        :return: A dictionary mapping each region to a list of oligos removed under this policy.
        :rtype: dict[str, list[str]]
        """
        oligos_with_hits = {region: [] for region in oligo_database.database.keys()}

        # remove all query oligos
        for hit in oligo_pair_hits:
            region_0 = self._oligo_to_region(hit[0])
            oligos_with_hits[str(region_0)].append(hit[0])

            region_1 = self._oligo_to_region(hit[1])
            oligos_with_hits[str(region_1)].append(hit[1])

        # remove duplicated entries
        for key in oligos_with_hits:
            oligos_with_hits[key] = list(set(oligos_with_hits[key]))

        return oligos_with_hits


class RemoveByLargerRegionPolicy(FilterPolicyBase):
    """A policy that removes oligonucleotides based on the number of oligos associated to a region.
    When two ologos (nodes in a graph) are connected by an edge, they are considered to potentially hybridize.
    For each edge, the oligo from the region with the higher number of associated oligos is removed,
    presuming a higher tolerance for loss in larger regions. This process is iterated until all edges are
    removed from the graph.
    """

    def __init__(self):
        """Constructor for the RemoveByLargerRegionPolicy class."""

    def apply(self, oligo_pair_hits: pd.DataFrame, oligo_database: OligoDatabase) -> dict[str, list[str]]:
        """Applies the policy to an oligo database, removing oligos to minimize
        hybridization based on the number of oligos associated with the region.

        :param oligo_pair_hits: Data frame containing pairs of oligonucleotides that potentially hybridize.
        :type oligo_pair_hits: pd.DataFrame
        :param oligo_database: The database of oligonucleotides.
        :type oligo_database: OligoDatabase
        :return: A dictionary mapping each region to a list of oligos removed under this policy.
        :rtype: dict[str, list[str]]
        """
        graph = nx.from_edgelist(oligo_pair_hits)
        number_oligos_per_region = self._get_number_oligos_per_region(oligo_database=oligo_database)
        oligos_with_hits = {region: [] for region in oligo_database.database.keys()}

        while graph.number_of_edges() > 0:
            edge = list(graph.edges)[0]
            region_0 = self._oligo_to_region(oligo=edge[0])
            region_1 = self._oligo_to_region(oligo=edge[1])
            if number_oligos_per_region[region_0] > number_oligos_per_region[region_1]:
                graph.remove_node(edge[0])
                oligos_with_hits[region_0].append(edge[0])
                number_oligos_per_region[region_0] -= 1
            else:
                graph.remove_node(edge[1])
                oligos_with_hits[region_1].append(edge[1])
                number_oligos_per_region[region_1] -= 1
        return oligos_with_hits


class RemoveByDegreePolicy(FilterPolicyBase):
    """A policy that removes oligonucleotides based on their node degree (number of connections to other nodes)
    in the graph. When two ologos (nodes in a graph) are connected by an edge, they are considered to potentially
    hybridize. Oligos with the highest node degree, i.e. the most interconnected nodes, are removed first.
    This process is iterated until all edges are removed from the graph.
    """

    def __init__(self):
        """Constructor for the RemoveByDegreePolicy class."""

    def apply(self, oligo_pair_hits: pd.DataFrame, oligo_database: OligoDatabase) -> dict[str, list[str]]:
        """Applies the policy to an oligo database, prioritizing the removal of oligos with the highest
        node degree of in the hybridization graph.

        :param oligo_pair_hits: Data frame containing pairs of oligonucleotides that potentially hybridize.
        :type oligo_pair_hits: pd.DataFrame
        :param oligo_database: The database of oligonucleotides.
        :type oligo_database: OligoDatabase
        :return: A dictionary mapping each region to a list of oligos removed under this policy.
        :rtype: dict[str, list[str]]
        """
        graph = nx.from_edgelist(oligo_pair_hits)
        oligos_with_hits = {region: [] for region in oligo_database.database.keys()}

        while graph.number_of_edges() > 0:
            degrees = dict(graph.degree())
            if degrees:
                max_degree_node = max(degrees, key=degrees.get)
                graph.remove_node(max_degree_node)
                oligos_with_hits[self._oligo_to_region(oligo=max_degree_node)].append(max_degree_node)
        return oligos_with_hits
