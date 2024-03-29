"""Extract each molecule of a GRO file and print their occurence.

Usage :
    import ModuleGroDecoder
"""

__authors__ = ("Karine DUONG")
__contact__ = ("karine.duong@etu.u-paris.fr")
__copyright__ = "IBPC"
__date__ = "2024-03-18"


import time
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis.distances import contact_matrix
from scipy.spatial.distance import cdist
import networkx as nx
from networkx.algorithms.components.connected import connected_components
import matplotlib.pyplot as plt
from loguru import logger


class Colors:
    """Text coloring code for terminal output."""

    RESET = '\033[0m'
    RED = '\033[31m'
    GREEN = '\033[32m'
    YELLOW = '\033[33m'
    BLUE = '\033[34m'
    MAGENTA = '\033[35m'
    CYAN = '\033[36m'
    WHITE = '\033[37m'
    BOLD = '\033[1m'

    def colorize_text(text, color):
        """Colorize text with the specified color."""
        return f"{color}{text}{Colors.RESET}"


BOND_LENGTH = {'C-C': 1.54,
               'C=C': 1.34,
               'C=O': 1.20,
               'O-H': 0.97,
               'C-H': 1.09,
               'N-H': 1.00,
               'C-S': 1.81,
               'S-H': 1.32,
               'N-C': 1.47,
               'C=N': 1.27,
               'S-S': 2.04,
               'C-O': 1.43
               }  # in Angstrom
# hydrogenBond = 2.7-3.3 A <--> 0.2-0.3 nm
# https://www.umass.edu/microbio/chime/find-ncb/help_gb.htm


def matrice(file_gro):
    """Calculate interatomic distances between all atoms in the GRO file \
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.cdist.html \
    https://stackoverflow.com/questions/72701992/convert-a-matrix-of-distance-to-adjacency-list/72702534#72702534 .

    (Warning : MDAnalysis convert the coordonate in Angstrom
    despite in the gro file it's in nm
    https://userguide.mdanalysis.org/stable/units.html#table-baseunits
    https://manual.gromacs.org/current/reference-manual/file-formats.html#gro )

    Args:
        file_gro (MDAnalysis.core.universe.Universe): An object
            representing the molecular structure loaded from a GRO file.

    Returns:
        numpy.ndarray : matrix of interatomic distances between all atoms
    """
    position = file_gro.atoms.positions
    taille = len(file_gro.atoms)
    tmp = np.full((taille, taille), 100.0)
    for i, pos_i in enumerate(position):
        for j, pos_j in enumerate(position[i+1:], start=i+1):
            tmp[i][j] = cdist([pos_i], [pos_j])
    return tmp


def pairmatrix_to_graph(liste):
    """Convert a list of pairs to a graph and its connected components.

    SOURCE : https://stackoverflow.com/questions/4842613/merge-lists-that-share-common-elements

    Args:
        liste (list): A list of pairs representing edges in the graph.

    Returns :
        tuple: A tuple containing two elements:
               1. A generator object that iterates over the connected components of the graph.
               2. A graph object representing the input list as a graph.
    """
    def to_graph(liste):
        """Create a graph based on a list of list \
           Where each sub-list is a group of connected node in the graph.

        Args:
            l (list): list of pair's list

        Returns:
            networkx.classes.graph.Graph: graph where each node is a value of the list
        """
        graph = nx.Graph()
        for part in liste:  # each sublist is a bunch of nodes
            graph.add_nodes_from(part)
            # it also imlies a number of edges:
            graph.add_edges_from(to_edges(part))
        return graph

    def to_edges(liste):
        """Generate edges for a graph from a list of nodes.

        Args :
            liste (list): A list of nodes representing a graph.
        """
        it = iter(liste)
        last = next(it)

        for current in it:
            yield last, current
            last = current

    graph = to_graph(liste)
    # networkx.draw(G, with_labels=True, font_weight='bold')
    # plt.show()
    return (connected_components(graph), graph)


def relabel_node(graph, G, atom_name):
    """Iterate each connected subgraph, and relabel their node by [atom_name]_[index].

    Args:
        graph (generator): An object that lists all connected graphs
        G (networkx.classes.graph.Graph): The graph for this file
        atom_name (numpy.ndarray): A list of each atoms for this file

    Returns:
        list: list where each node is relabel
    """
    modified_subgraphs = []   # List to store modified subgraphs

    # For each connected component, create a subgraph and add it to the list
    for component in graph:
        subgraph = G.subgraph(component)  # Create a subgraph from the component

        # Create a replacement dictionary for node names
        mapping = {node: f"{atom_name[node]}_{node}" for node in subgraph.nodes()}
        # Create a new subgraph with modified node names
        modified_subgraph = nx.relabel_nodes(subgraph, mapping)
        modified_subgraphs.append(modified_subgraph)  # Add the modified subgraph to the list of modified subgraphs
    return modified_subgraphs


def extract_node_labels(graph):
    """Extract node labels (just the atom's name) from a graph.

    Args:
        graph (networkx.Graph): The input graph.

    Returns:
        tuple: A tuple containing the node labels extracted from the graph sorted.
    """
    return tuple(sorted([str(node).split('_')[0] for node in graph.nodes()]))


def compare_node_labels(graph1, graph2):
    """Compare the node labels of two graphs.

    Args:
        graph1 (networkx.Graph): The first input graph.
        graph2 (networkx.Graph): The second input graph.

    Returns:
        bool: True if the node labels of both graphs are the same, False otherwise.
    """
    labels1 = extract_node_labels(graph1)
    labels2 = extract_node_labels(graph2) 
    return labels1 == labels2


def compare_graphs(graph1, graph2):
    """Compare two graphs for isomorphism.

    Args:
        graph1 (networkx.Graph): The first input graph.
        graph2 (networkx.Graph): The second input graph.

    Returns:
        bool: True if the graphs are isomorphic, False otherwise.
    """
    if not compare_node_labels(graph1, graph2):  # first see if both graph have the same label
        return False
    else:  # if it's True, see if both graph have the same connection
        return nx.vf2pp_is_isomorphic(graph1, graph2) #node_label="label"


def count_molecule(modified_subgraph):
    """Print information about the modified subgraphs.

    Args:
        modified_subgraph (list): A list of modified subgraphs.

    Returns:
        dict: A dictionary containing information about the subgraphs.
              Keys are tuples representing the node labels of each subgraph,
              and values are the number of occurrences of each subgraph.
    """
    size = {}  # Dictionary to store each subgraph based on their size
    count = {}  # Dictionary to store each subgraph based on their occurence
    for subgraph in modified_subgraph:  # Iterate on each sublist, and sort it based on their size
        taille = subgraph.size()
        if taille in size:
            size[taille].append(subgraph)
        else:
            size[taille] = [subgraph]

    for taille, subgraph_list in size.items():  # Iterate on the dictionary
        k = 0
        while k < len(subgraph_list):
            if len(subgraph_list) == 1:  # If the number of subgraph for this size is 1
                name_molecule = tuple(extract_node_labels(subgraph_list[k]))
                count[name_molecule] = 1  # Add it to the dictionary
            else : 
                if taille==1 :  # See if the subgraph only have one node
                    name_node = next(iter(subgraph_list[k].nodes()))
                    nom = str(name_node).split('_')[0]
                    count[nom] = count.get(nom, 0) + 1 
                else : 
                    i = k 
                    while i < len(subgraph_list)-1 : 
                        if compare_graphs(subgraph_list[i], subgraph_list[i+1]): 
                            name_node = extract_node_labels(subgraph_list[i])
                            count[name_node] = count.get(name_node, 1) + 1 
                        else : 
                            name_node1 = extract_node_labels(subgraph_list[i])
                            name_node2 = extract_node_labels(subgraph_list[i+1])
                            count[name_node1] = count.get(name_node1, 0) + 1 
                            count[name_node2] = count.get(name_node2, 0) + 1 
                        i += 1 
                        k = i+1
            k += 1
    return count


def print_count(count, option=""):
    """Print the count of molecules.

    Args:
        count (dict): A dictionary containing information about the molecules.
        option (str): either we want to print the composition of each molecule (with the option "detail") or not
    """
    logger.info("Voici le contenu de ce fichier GRO : ")
    size = 0
    for i, (key, value) in enumerate(count.items()):
        logger.info(f"\nMolecule {i+1}: \n\t {len(key):,} atoms \n\t Quantity: {value:,}")
        if option=="detail":
            logger.info(f"\t Composition: {key}")
        size += value
    logger.info(f"It containt {size:,} molecule")


def name_molecule_graph(belong_mole, file_gro):
    """Based on the value of each list (one molecule), it search the name of each atome.

    Args:
        belong_mole (list): list of every indice (from de matrix) that have bond between them
        file_gro (MDAnalysis.core.universe.Universe): object that contain information about GRO file

    Returns:
        list: list that contain the name of each atom for each molecule (sub-list)
    """
    atoms_names = file_gro.atoms.names
    list_mole = []

    for sublist in belong_mole:
        mol = [atoms_names[indice] for indice in sublist]
        list_mole.append(mol)
    return list_mole


def delete_hydrogen_grofile(filepath):
    """Delete hydrogen atoms from a .gro file and save the modified file.

    Args:
        filepath (str): The filepath of the .gro file containing hydrogen atoms.

    Returns:
        Universe: A modified MDA Universe object containing only non-hydrogen atoms.
    """
    filegro = mda.Universe(filepath)
    atoms_remains = filegro.select_atoms("not (name H*)")
    new_filepath = f"{filepath[:-4]}_woH.gro"
    atoms_remains.write(new_filepath)
    return atoms_remains


def grodecoder_principal(filepath_gro):
    """Excute the main function for analyzing a .gro file.

    Args:
        filepath_gro (str): Filepath of the .gro file we want to analyzed
    """
    logger.info(f"Filename: {filepath_gro} --------") 

    threshold = max(BOND_LENGTH.values())
    logger.info(f"Threshold: {threshold} Angstrom")

    file_gro = mda.Universe(filepath_gro)  # load the .gro file
    logger.info(f"How many atoms there is in this file (at the beginning): {len(file_gro.atoms):,}")

    file_gro = delete_hydrogen_grofile(filepath_gro)
    logger.info(f"How many atoms there is in this file (without all hydrogen): {len(file_gro.atoms):,}")

    # https://docs.mdanalysis.org/1.1.0/documentation_pages/analysis/distances.html#MDAnalysis.analysis.distances.contact_matrix
    start_time = time.time()
    mat_contact = contact_matrix(file_gro.atoms.positions, cutoff=threshold, returntype='sparse')
    end_time = time.time()
    temps = end_time - start_time
    logger.success(Colors.colorize_text(f"[contact_matrix]{temps} seconds", Colors.BLUE))

    # https://numpy.org/doc/stable/reference/generated/numpy.argwhere.html
    # https://www.includehelp.com/python/how-to-get-indices-of-elements-that-are-greater-than-a-threshold-in-2d-numpy-array.aspx
    start_time = time.time()
    pair_matrix = np.argwhere(mat_contact)  # list of pair where their value (their indice in the distance matrix) is below the threshold
    end_time = time.time()
    temps = end_time - start_time
    logger.success(Colors.colorize_text(f"[pair_matrix]{temps} seconds", Colors.BLUE))

    # Turn each tuple into a graph, if one node's label is already existant
    # It attached the other label to it
    # So it connected all the node who have common label
    start_time = time.time()
    connexgraph_return, graph_return = pairmatrix_to_graph(pair_matrix)
    modified_subgraph = relabel_node(connexgraph_return, graph_return, file_gro.atoms.names)
    count = count_molecule(modified_subgraph)
    print_count(count, "detail")
    end_time = time.time()
    temps = end_time - start_time
    logger.success(Colors.colorize_text(f"[pairmatrix_to_graph / count_molecule / print_graph / print_count ]{temps} seconds", Colors.BLUE))
    print("[main] ------------")

