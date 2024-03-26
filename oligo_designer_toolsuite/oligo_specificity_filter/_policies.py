############################################
# imports
############################################

from abc import ABC, abstractmethod

import networkx as nx

from .._constants import SEPARATOR_OLIGO_ID
from ..database import OligoDatabase

############################################
# Crosshybridization Policies
############################################


class FilterPolicyBase(ABC):
    """An abstract base class for defining the policy for managing cross-hybridization within an oligonucleotide database.
    Implementations of this class should provide specific strategies to minimize cross-hybridization risk in
    the design of oligonucleotide sequences.
    """

    def __init__(self):
        """Constructor for the CrossHybridizationPolicy class."""

    @abstractmethod
    def apply(self, graph: nx.Graph, oligo_database: OligoDatabase):
        """Applies the cross-hybridization policy to an oligonucleotide database, removing oligos to reduce
        cross-hybridization risks.

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


class RemoveByLargerRegionPolicy(FilterPolicyBase):
    """A policy that removes oligonucleotides based on the number of oligos associated to a regions to reduce cross-hybridization.
    When two ologos (nodes in a graph) are connected by an edge, they are considered to potentially cross-hybridize.
    For each edge, the oligo from the region with the higher number of associated oligos is removed, presuming a higher tolerance for loss in larger regions.
    This process is iterated until all edges are removed from the graph.
    """

    def __init__(self):
        """Constructor for the RemoveByLargerRegionPolicy class."""

    def apply(self, graph: nx.Graph, oligo_database: OligoDatabase) -> dict[str, list[str]]:
        """Applies the policy to an oligo database, removing oligos to minimize cross-hybridization based on the number of oligos associated with the region.

        :param graph: A graph where nodes represent oligos and edges indicate potential cross-hybridization.
        :type graph: nx.Graph
        :param oligo_database: The database of oligonucleotides.
        :type oligo_database: OligoDatabase
        :return: A dictionary mapping each region to a list of oligos removed under this policy.
        :rtype: dict[str, list[str]]
        """
        number_oligos_per_region = self._get_number_oligos_per_region(oligo_database=oligo_database)
        removed_oligos = {region: [] for region in oligo_database.database.keys()}

        while graph.number_of_edges() > 0:
            edge = list(graph.edges)[0]
            region_0 = self._oligo_to_region(oligo=edge[0])
            region_1 = self._oligo_to_region(oligo=edge[1])
            if number_oligos_per_region[region_0] > number_oligos_per_region[region_1]:
                graph.remove_node(edge[0])
                removed_oligos[region_0].append(edge[0])
                number_oligos_per_region[region_0] -= 1
            else:
                graph.remove_node(edge[1])
                removed_oligos[region_1].append(edge[1])
                number_oligos_per_region[region_1] -= 1
        return removed_oligos


class RemoveByDegreePolicy(FilterPolicyBase):
    """A policy that removes oligonucleotides based on their node degree (number of connections to other nodes) in the cross-hybridization graph.
    When two ologos (nodes in a graph) are connected by an edge, they are considered to potentially cross-hybridize.
    Oligos with the highest node degree, i.e. the most interconnected nodes, are removed first. This process is iterated until all edges are removed from the graph.
    """

    def __init__(self):
        """Constructor for the RemoveByDegreePolicy class."""

    def apply(self, graph: nx.Graph, oligo_database: OligoDatabase) -> dict[str, list[str]]:
        """Applies the policy to an oligo database, prioritizing the removal of oligos with the highest node degree of in the cross-hybridization graph.

        :param graph: A graph where nodes represent oligos and edges indicate potential cross-hybridization.
        :type graph: nx.Graph
        :param oligo_database: The database of oligonucleotides.
        :type oligo_database: OligoDatabase
        :return: A dictionary mapping each region to a list of oligos removed under this policy.
        :rtype: dict[str, list[str]]
        """
        removed_oligos = {region: [] for region in oligo_database.database.keys()}

        while graph.number_of_edges() > 0:
            degrees = dict(graph.degree())
            if degrees:
                max_degree_node = max(degrees, key=degrees.get)
                graph.remove_node(max_degree_node)
                removed_oligos[self._oligo_to_region(oligo=max_degree_node)].append(max_degree_node)
        return removed_oligos
