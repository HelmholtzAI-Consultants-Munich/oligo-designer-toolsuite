############################################
# imports
############################################

import os
from abc import ABC, abstractmethod
from pathlib import Path

import networkx as nx

from .._constants import _TYPES_SEQ, SEPARATOR_OLIGO_ID
from ..database import OligoDatabase, ReferenceDatabase
from . import AlignmentSpecificityFilter, SpecificityFilterBase

############################################
# Crosshybridization Policies
############################################


class CrossHybridizationPolicy(ABC):
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


class RemoveByLargerRegionPolicy(CrossHybridizationPolicy):
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


class RemoveByDegreePolicy(CrossHybridizationPolicy):
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


############################################
# Oligo Crosshybridization Filter Classe
############################################


class CrossHybridizationFilter(SpecificityFilterBase):
    """A class that implements a filter to minimize cross-hybridization within an oligonucleotide database by applying a specified policy.
    This class integrates a cross-hybridization policy and a specificity filter to identify and mitigate potential
    cross-hybridization events between oligonucleotides.

    :param policy: The cross-hybridization policy to apply for minimizing cross-hybridization.
    :type policy: CrossHybridizationPolicy
    :param specificity_filter: The alignment specificity filter used to identify potential cross-hybridization events.
    :type specificity_filter: AlignmentSpecificityFilter
    :param dir_cross_hybridization: Directory for saving output files related to cross-hybridization filtering.
    :type dir_cross_hybridization: str
    """

    def __init__(
        self,
        policy: CrossHybridizationPolicy,
        specificity_filter: AlignmentSpecificityFilter,
        dir_output: str = "output",
    ):
        """Constructor for the CrossHybridizationFilter class."""
        self.policy = policy
        self.specificity_filter = specificity_filter

        self.dir_cross_hybridization = os.path.join(dir_output, "crosshybridization")
        Path(self.dir_cross_hybridization).mkdir(parents=True, exist_ok=True)

    def apply(self, sequence_type: _TYPES_SEQ, oligo_database: OligoDatabase, n_jobs: int):
        """Applies the cross-hybridization filter to an oligonucleotide database.

        :param sequence_type: The type of sequences being filtered, must be one of the predefined sequence types.
        :type sequence_type: _TYPES_SEQ
        :param oligo_database: The database of oligonucleotides to be filtered.
        :type oligo_database: OligoDatabase
        :param n_jobs: The number of parallel jobs to run.
        :type n_jobs: int
        :return: The oligo database with cross-hybridization minimized according to the policy.
        :rtype: OligoDatabase
        """
        regions = list(oligo_database.database.keys())

        reference_database = self._create_reference_database(
            sequence_type=sequence_type, oligo_database=oligo_database
        )
        cross_hybridization_graph = self._create_cross_hybridization_graph(
            sequence_type=sequence_type,
            oligo_database=oligo_database,
            reference_database=reference_database,
            n_jobs=n_jobs,
        )
        oligos_with_hits = self.policy.apply(graph=cross_hybridization_graph, oligo_database=oligo_database)

        for region in regions:
            database_region_filtered = self._filter_hits_from_database(
                database_region=oligo_database.database[region],
                oligos_with_hits=oligos_with_hits[region],
            )
            oligo_database.database[region] = database_region_filtered
        return oligo_database

    def _create_reference_database(self, sequence_type: _TYPES_SEQ, oligo_database: OligoDatabase):
        """Creates a reference database from the given oligo database. This involves writing the oligo database to a FASTA file and loading it into a new ReferenceDatabase instance, including metadata.

        :param sequence_type: The type of sequences being filtered, must be one of the predefined sequence types.
        :type sequence_type: _TYPES_SEQ
        :param oligo_database: The oligo database from which the reference database is created.
        :type oligo_database: OligoDatabase
        :return: A ReferenceDatabase instance populated with sequences and metadata from the oligo database.
        :rtype: ReferenceDatabase
        """
        file_reference = oligo_database.write_database_to_fasta(
            filename=f"oligo_database_crosshybridization_with_{sequence_type}",
            region_ids=None,
            sequence_type=sequence_type,
        )
        reference_database = ReferenceDatabase()
        reference_database.load_metadata(metadata=oligo_database.metadata)
        reference_database.load_sequences_fom_fasta(file_fasta=file_reference, database_overwrite=True)

        os.remove(file_reference)

        return reference_database

    def _create_cross_hybridization_graph(
        self,
        sequence_type: _TYPES_SEQ,
        oligo_database: OligoDatabase,
        reference_database: ReferenceDatabase,
        n_jobs: int,
    ):
        """Creates a graph representing potential cross-hybridization events between oligonucleotides.

        :param sequence_type: The type of sequences being filtered, must be one of the predefined sequence types.
        :type sequence_type: _TYPES_SEQ
        :param oligo_database: The oligonucleotide database.
        :type oligo_database: OligoDatabase
        :param reference_database: The reference database to compare against for cross-hybridization.
        :type reference_database: ReferenceDatabase
        :param n_jobs: The number of parallel jobs to run.
        :type n_jobs: int
        :return: A graph where nodes represent oligonucleotides and edges represent potential cross-hybridization.
        :rtype: nx.Graph
        """
        oligo_pairs_hits = self.specificity_filter.get_oligo_pair_hits(
            sequence_type=sequence_type,
            oligo_database=oligo_database,
            n_jobs=n_jobs,
            reference_database=reference_database,
        )

        return nx.from_edgelist(oligo_pairs_hits)
