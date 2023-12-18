import networkx as nx

from oligo_designer_toolsuite.constants import REGION_SEPARATOR_OLIGO


# Policies
def _oligo_to_region(oligo_name: str):
    return oligo_name.split(REGION_SEPARATOR_OLIGO)[0]


def remove_by_bigger_region_policy(
    graph: nx.Graph, oligos_per_region: dict
) -> dict[str, list[str]]:
    """
    Iteratively removes oligos (graph nodes) by comparing the number of oligos per region. In each iteration, the oligo from the region with the bigger number of oligos is removed. The process is repeated until there are no edges in the graph anymore.

    :param graph: Graph where nodes represent oligos and edges indicate cross-hybridization between pairs of oligos.
    :type graph: nx.Graph
    :param oligos_per_region: Dictionary mapping each region to its number of oligos. Nodes from regions with more oligos are prioritized for removal.
    :type oligos_per_region: dict
    :return: Dictionary where each key is a region and the value is a list of oligos removed from that region.
    :rtype: dict
    """

    removed_oligos = {oligo_region: [] for oligo_region in oligos_per_region.keys()}

    while graph.number_of_edges() > 0:
        edge = list(graph.edges)[0]
        region_0 = _oligo_to_region(edge[0])
        region_1 = _oligo_to_region(edge[1])
        if oligos_per_region[region_0] > oligos_per_region[region_1]:
            graph.remove_node(edge[0])
            removed_oligos[region_0].append(edge[0])
            oligos_per_region[region_0] -= 1
        else:
            graph.remove_node(edge[1])
            removed_oligos[region_1].append(edge[1])
            oligos_per_region[region_1] -= 1
    return removed_oligos


def remove_by_degree_policy(
    graph: nx.Graph, oligos_per_region: dict
) -> dict[str, list[str]]:
    """
    Iteratively removes nodes with the highest degree from a graph until there are no edges left. The degree of a node is the number of connections it has with other nodes. This method targets nodes (oligos) that are most interconnected first.

    :param graph: Graph where nodes represent oligos and edges indicate cross-hybridization between pairs of oligos.
    :type graph: nx.Graph
    :param oligos_per_region: Dictionary mapping each region to its number of oligos.
    :type oligos_per_region: dict
    :return: Dictionary where each key is a region and the value is a list of oligos removed from that region.
    :rtype: dict
    """
    removed_oligos = {oligo_region: [] for oligo_region in oligos_per_region.keys()}

    while graph.number_of_edges() > 0:
        degrees = dict(graph.degree())
        if degrees:
            max_degree_node = max(degrees, key=degrees.get)
            graph.remove_node(max_degree_node)
            removed_oligos[_oligo_to_region(max_degree_node)].append(max_degree_node)
    return removed_oligos
