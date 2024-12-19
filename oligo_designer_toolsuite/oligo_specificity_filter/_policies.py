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
    def __init__(self) -> None:
        """Constructor for the FilterPolicyBase class."""

    @abstractmethod
    def apply(self, oligo_pair_hits: pd.DataFrame, oligo_database: OligoDatabase) -> dict:
        """ """

    def _oligo_to_region(self, oligo: str) -> str:
        """ """
        return oligo.split(SEPARATOR_OLIGO_ID)[0]

    def _get_number_oligos_per_region(self, oligo_database: OligoDatabase) -> dict:
        """ """
        return {region: len(oligo_database.database[region]) for region in oligo_database.database.keys()}


class RemoveAllPolicy(FilterPolicyBase):
    def __init__(self) -> None:
        """Constructor for the RemoveAllPolicy class."""

    def apply(self, oligo_pair_hits: pd.DataFrame, oligo_database: OligoDatabase) -> dict:
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
    def __init__(self) -> None:
        """Constructor for the RemoveByLargerRegionPolicy class."""

    def apply(self, oligo_pair_hits: pd.DataFrame, oligo_database: OligoDatabase) -> dict:
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
    def __init__(self) -> None:
        """Constructor for the RemoveByDegreePolicy class."""

    def apply(self, oligo_pair_hits: pd.DataFrame, oligo_database: OligoDatabase) -> dict:
        graph = nx.from_edgelist(oligo_pair_hits)
        oligos_with_hits = {region: [] for region in oligo_database.database.keys()}

        while graph.number_of_edges() > 0:
            degrees = dict(graph.degree())
            if degrees:
                max_degree_node = max(degrees, key=degrees.get)
                graph.remove_node(max_degree_node)
                oligos_with_hits[self._oligo_to_region(oligo=max_degree_node)].append(max_degree_node)
        return oligos_with_hits
