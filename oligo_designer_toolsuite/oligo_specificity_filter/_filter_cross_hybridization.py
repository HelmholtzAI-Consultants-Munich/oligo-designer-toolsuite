############################################
# imports
############################################

import networkx as nx

from . import AlignmentSpecificityFilter, SpecificityFilterBase

############################################
# Oligo Blast Filter Classes
############################################


class CrossHybridizationFilter(SpecificityFilterBase):
    """
    This class implements a cross-hybridization filter for oligos. It constructs a graph representation of oligo interactions using an alignment method, then applies a policy to remove oligos which cross-hybridize.
    The class provides methods to create a cross-hybridization graph, and apply the chosen policy to filter out oligos.
    The policy is a function that determines which oligos to remove in case of cross-hybridization (can be imported from oligo_designer_toolsuite.oligo_specificity_filter.cross_hybridization_policies, or you can provide your own).

    :param specificity_filter: An instance of AlignmentSpecificityFilter class that provides matching oligo pairs based on an alignment criteria.
    :type specificity_filter: AlignmentSpecificityFilter
    :param policy: A function that defines the strategy for oligo removal. It takes the graph and a dictionary where keys are the databse regions as inputs and returns a dictionary of removed oligos.
    :type policy: function
    :param dir_cross_hybridization: Directory where files produced during the cross hybridization check are stored.
    :type dir_cross_hybridization: str
    """

    def __init__(
        self,
        specificity_filter: AlignmentSpecificityFilter,
        policy,
        dir_cross_hybridization,
    ):
        self.specificity_filter = specificity_filter
        self.policy = policy
        self.dir_cross_hybridization = dir_cross_hybridization

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
        oligos_per_region = self._count_oligos_per_region(database)
        regions = list(database.keys())

        cross_hybridization_graph = self._create_cross_hybridization_graph(
            database, file_reference
        )
        matching_oligos = self.policy(cross_hybridization_graph, oligos_per_region)

        for region in regions:
            filtered_database_region = self._filter_matching_oligos(
                database[region],
                matching_oligos[region],
            )
            database[region] = filtered_database_region
        return database

    def _create_cross_hybridization_graph(
        self, oligo_database: dict, reference_fasta: str
    ):
        matching_oligo_pairs = self.specificity_filter.get_matching_oligo_pairs(
            oligo_database, reference_fasta
        )
        print(matching_oligo_pairs)
        return nx.from_edgelist(matching_oligo_pairs)

    def _count_oligos_per_region(self, oligo_database):
        return {region: len(oligo_database[region]) for region in oligo_database.keys()}
