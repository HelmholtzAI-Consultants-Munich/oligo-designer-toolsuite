from abc import ABC, abstractmethod

import networkx as nx

from oligo_designer_toolsuite.constants import REGION_SEPARATOR_OLIGO
from oligo_designer_toolsuite.database import OligoDatabase


class CrossHybridizationPolicy(ABC):
    """
    This is the base class for all cross-hybridization policies.
    """

    @abstractmethod
    def apply(self, graph: nx.Graph, oligo_database: OligoDatabase):
        """
        Apply the cross-hybridization policy to a graph of oligos.
        :param graph: Graph where nodes represent oligos and edges indicate cross-hybridization between pairs of oligos.
        :type graph: nx.Graph
        :param oligo_database: The oligo database.
        :type oligo_database: OligoDatabase

        :return: Dictionary where each key is a region and the value is a list of oligos removed from that region.
        :rtype: dict
        """

    def _oligo_to_region(self, oligo: str) -> str:
        """
        Extract the region from an oligo ID.
        :param oligo: The oligo ID.
        :type oligo: str
        :return: The region.
        :rtype: str
        """
        return oligo.split(REGION_SEPARATOR_OLIGO)[0]

    def _count_oligos_per_region(self, oligo_database):
        """
        Count the number of oligos per region in the provided oligo database.

        :param oligo_database: A dictionary containing oligos grouped by regions.
        :type oligo_database: dict
        :return: A dictionary mapping region names to the number of oligos in each region.
        :rtype: dict
        """
        return {region: len(oligo_database[region]) for region in oligo_database.keys()}


class RemoveByBiggerRegionPolicy(CrossHybridizationPolicy):
    """
    Iteratively removes oligos (graph nodes) by comparing the number of oligos per region. In each iteration, the oligo from the region with the bigger number of oligos is removed. The process is repeated until there are no edges in the graph anymore.
    """

    def apply(
        self, graph: nx.Graph, oligo_database: OligoDatabase
    ) -> dict[str, list[str]]:
        """
        Apply the cross-hybridization policy to a graph of oligos.
        :param graph: Graph where nodes represent oligos and edges indicate cross-hybridization between pairs of oligos.
        :type graph: nx.Graph
        :param oligo_database: The oligo database.
        :type oligo_database: OligoDatabase

        :return: Dictionary where each key is a region and the value is a list of oligos removed from that region.
        :rtype: dict
        """
        oligos_per_region = self._count_oligos_per_region(oligo_database)
        removed_oligos = {oligo_region: [] for oligo_region in oligos_per_region.keys()}

        while graph.number_of_edges() > 0:
            edge = list(graph.edges)[0]
            region_0 = self._oligo_to_region(edge[0])
            region_1 = self._oligo_to_region(edge[1])
            if oligos_per_region[region_0] > oligos_per_region[region_1]:
                graph.remove_node(edge[0])
                removed_oligos[region_0].append(edge[0])
                oligos_per_region[region_0] -= 1
            else:
                graph.remove_node(edge[1])
                removed_oligos[region_1].append(edge[1])
                oligos_per_region[region_1] -= 1
        return removed_oligos


class RemoveByDegreePolicy(CrossHybridizationPolicy):
    """
    Iteratively removes nodes with the highest degree from a graph until there are no edges left. The degree of a node is the number of connections it has with other nodes. This method targets nodes (oligos) that are most interconnected first.
    """

    def apply(
        self, graph: nx.Graph, oligo_database: OligoDatabase
    ) -> dict[str, list[str]]:
        """
        Apply the cross-hybridization policy to a graph of oligos.
        :param graph: Graph where nodes represent oligos and edges indicate cross-hybridization between pairs of oligos.
        :type graph: nx.Graph
        :param oligo_database: The oligo database.
        :type oligo_database: OligoDatabase

        :return: Dictionary where each key is a region and the value is a list of oligos removed from that region.
        :rtype: dict
        """
        removed_oligos = {oligo_region: [] for oligo_region in oligo_database.keys()}

        while graph.number_of_edges() > 0:
            degrees = dict(graph.degree())
            if degrees:
                max_degree_node = max(degrees, key=degrees.get)
                graph.remove_node(max_degree_node)
                removed_oligos[self._oligo_to_region(max_degree_node)].append(
                    max_degree_node
                )
        return removed_oligos
