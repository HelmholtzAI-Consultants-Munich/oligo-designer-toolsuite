############################################
# imports
############################################

import os

import networkx as nx

from . import AlignmentSpecificityFilter, SpecificityFilterBase
from .cross_hybridization_policies import CrossHybridizationPolicy

############################################
# Oligo Blast Filter Classes
############################################


class CrossHybridizationFilter(SpecificityFilterBase):
    """
    This class implements a cross-hybridization filter for oligos. It constructs a graph representation of oligo interactions using an alignment method, then applies a policy to remove oligos which cross-hybridize.
    The class provides methods to create a cross-hybridization graph, and apply the chosen policy to filter out oligos.
    The policy determines which oligos to remove in case of cross-hybridization.

    :param specificity_filter: An instance of AlignmentSpecificityFilter class that provides matching oligo pairs based on an alignment criteria.
    :type specificity_filter: AlignmentSpecificityFilter
    :param policy: Strategy for oligo removal. Can be imported from oligo_designer_toolsuite.oligo_specificity_filter.cross_hybridization_policies, or you can provide your own.
    :type policy: CrossHybridizationPolicy
    :param dir_cross_hybridization: Directory where files produced during the cross hybridization check are stored.
    :type dir_cross_hybridization: str
    """

    def __init__(
        self,
        specificity_filter: AlignmentSpecificityFilter,
        policy: CrossHybridizationPolicy,
        dir_cross_hybridization: str,
    ):
        self.specificity_filter = specificity_filter
        self.policy = policy
        if dir_cross_hybridization is not None:
            self.dir_cross_hybridization = dir_cross_hybridization
        else:
            self.dir_cross_hybridization = os.path.join(
                specificity_filter.dir_specificity, "cross_hybridization"
            )

    def apply(self, database: dict, file_reference: str, n_jobs: int):
        """
        Applies the cross-hybridization filter to an oligo database. First, it constructs a graph of oligo interactions. Then, it identifies cross-hybridizing oligos, and filters one of them out based on the specified policy.

        :param database: The oligo database to be filtered.
        :type database: dict
        :param file_reference: Path to the reference file containing all oligos in fasta.
        :type file_reference: str
        :param n_jobs: Number of parallel jobs to use for parallel processing (currently not implemented).
        :type n_jobs: int
        :return: Updated oligo database with cross-hybridizing oligos removed.
        :rtype: dict
        """
        # Not in parallel for now
        regions = list(database.keys())

        cross_hybridization_graph = self._create_cross_hybridization_graph(
            database, file_reference, n_jobs
        )
        matching_oligos = self.policy.apply(cross_hybridization_graph, database)

        for region in regions:
            filtered_database_region = self._filter_matching_oligos(
                database[region],
                matching_oligos[region],
            )
            database[region] = filtered_database_region
        return database

    def _create_cross_hybridization_graph(
        self, oligo_database: dict, reference_fasta: str, n_jobs: int
    ):
        """
        Create a cross-hybridization graph based on matching oligo pairs. Where nodes are oligos and edges indicate that the oligos match

        :param oligo_database: The oligo database
        :type oligo_database: dict
        :param reference_fasta: The path to the reference FASTA file.
        :type reference_fasta: str
        :param n_jobs: Number of parallel jobs to use for parallel processing.
        :type n_jobs: int
        :return: A graph representing cross-hybridization relationships between oligos.
        :rtype: networkx.Graph
        """
        matching_oligo_pairs = self.specificity_filter.get_matching_oligo_pairs(
            oligo_database, reference_fasta, n_jobs
        )
        print(matching_oligo_pairs)
        return nx.from_edgelist(matching_oligo_pairs)
