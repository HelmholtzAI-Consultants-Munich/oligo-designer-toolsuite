############################################
# imports
############################################

from abc import ABC, abstractmethod

import networkx as nx
import pandas as pd

from oligo_designer_toolsuite._constants import SEPARATOR_OLIGO_ID
from oligo_designer_toolsuite.database import OligoDatabase

############################################
# Hybridization Policies
############################################


class FilterPolicyBase(ABC):
    """
    An abstract base class for defining filter policies used in specficity filters.

    The `FilterPolicyBase` class provides the foundational structure for creating custom filtering policies
    that dictate how oligonucleotide hits are processed and filtered in an OligoDatabase.
    It includes essential methods for applying the filter and utility functions for managing and organizing oligos by region.

    """

    def __init__(self) -> None:
        """Constructor for the FilterPolicyBase class."""

    @abstractmethod
    def apply(self, oligo_pair_hits: pd.DataFrame, oligo_database: OligoDatabase) -> dict:
        """
        Abstract method to apply the filter policy on the oligo pair hits.

        :param oligo_pair_hits: DataFrame containing pairs of oligonucleotides that have been identified as hits.
        :type oligo_pair_hits: pd.DataFrame
        :param oligo_database: The OligoDatabase containing the oligonucleotides and their associated attributes.
        :type oligo_database: OligoDatabase
        :return: A dictionary mapping each region to a list of oligos that should be removed based on the policy.
        :rtype: dict
        """

    def _oligo_to_region(self, oligo: str) -> str:
        """
        Helper method to extract the region ID from an oligo ID.

        :param oligo: The oligo ID string from which the region ID will be extracted.
        :type oligo: str
        :return: The extracted region ID.
        :rtype: str
        """
        return oligo.split(SEPARATOR_OLIGO_ID)[0]

    def _get_number_oligos_per_region(self, oligo_database: OligoDatabase) -> dict:
        """
        Helper method to get the number of oligos in each region of the oligo database.

        :param oligo_database: The OligoDatabase containing the oligonucleotides and their associated attributes.
        :type oligo_database: OligoDatabase
        :return: A dictionary mapping each region to the number of oligos it contains.
        :rtype: dict
        """
        return {region: len(oligo_database.database[region]) for region in oligo_database.database.keys()}


class RemoveAllPolicy(FilterPolicyBase):
    """
    A filter policy that removes all oligonucleotides involved in any hits.

    The `RemoveAllPolicy` class identifies all oligonucleotides involved in a hit pair and removes them from the database.
    This policy ensures that any oligos associated with potential hybridization or other conflicts are excluded from further analysis.
    """

    def __init__(self) -> None:
        """Constructor for the RemoveAllPolicy class."""

    def apply(self, oligo_pair_hits: pd.DataFrame, oligo_database: OligoDatabase) -> dict:
        """
        Applies the filter policy by identifying all oligonucleotides involved in hits and marking them for removal.

        The `apply` method processes each pair of oligonucleotide hits and associates them with their respective regions.
        It then compiles a list of oligos to be removed based on these hits.

        :param oligo_pair_hits: DataFrame containing pairs of oligonucleotides that have been identified as hits.
        :type oligo_pair_hits: pd.DataFrame
        :param oligo_database: The OligoDatabase containing the oligonucleotides and their associated attributes.
        :type oligo_database: OligoDatabase
        :return: A dictionary mapping each region to a list of oligos that should be removed based on the policy.
        :rtype: dict
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
    """
    A filter policy class that removes oligonucleotides based on the number of oligos attributed to the region they belong to.

    The `RemoveByLargerRegionPolicy` class applies a filtering strategy where oligonucleotides involved in cross-region hits are
    removed based on the number of oligos in their respective regions. The oligo from the larger region (with more oligos) is removed first.
    """

    def __init__(self) -> None:
        """Constructor for the RemoveByLargerRegionPolicy class."""

    def apply(self, oligo_pair_hits: pd.DataFrame, oligo_database: OligoDatabase) -> dict:
        """
        Applies the policy to remove oligonucleotides based on number of oligos per region.

        This function processes pairs of oligonucleotide hits, building a graph from the edges representing these pairs.
        It iteratively removes oligos from larger regions (regions with more oligos) until no more edges remain in the graph.
        The removed oligos are tracked and returned in a dictionary keyed by region.

        :param oligo_pair_hits: DataFrame containing pairs of oligonucleotides that have been identified as hits.
        :type oligo_pair_hits: pd.DataFrame
        :param oligo_database: The OligoDatabase containing the oligonucleotides and their associated attributes.
        :type oligo_database: OligoDatabase
        :return: A dictionary mapping each region to a list of oligos that should be removed based on the policy.
        :rtype: dict
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
    """
    A filtering policy that removes oligonucleotides based on their connectivity within a graph of oligo pair hits.

    The `RemoveByDegreePolicy` class is designed to eliminate oligonucleotides that are most connected (i.e., those with the highest degree)
    within a graph of oligo pair hits. This approach aims to reduce the potential for hybridization by iteratively removing the most problematic
    sequences until no hybridization edges remain in the graph.
    """

    def __init__(self) -> None:
        """Constructor for the RemoveByDegreePolicy class."""

    def apply(self, oligo_pair_hits: pd.DataFrame, oligo_database: OligoDatabase) -> dict:
        """
        Applies the degree-based filtering policy by removing oligonucleotides with the highest degree of connectivity.

        :param oligo_pair_hits: DataFrame containing pairs of oligonucleotides that have been identified as hits.
        :type oligo_pair_hits: pd.DataFrame
        :param oligo_database: The OligoDatabase containing the oligonucleotides and their associated attributes.
        :type oligo_database: OligoDatabase
        :return: A dictionary mapping each region to a list of oligos that should be removed based on the policy.
        :rtype: dict
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
