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
from itertools import groupby 


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

def convert_atompairs_to_graph(atom_pairs, number_of_atoms):
    """Convert a list of pairs to a graph and its connected components.

    SOURCE : https://stackoverflow.com/questions/4842613/merge-lists-that-share-common-elements

    Parameters
    ----------
        atom_pairs : list
            A list of pairs representing edges in the graph.
        number_of_atoms : 
            

    Returns
    -------
        tuple
            A tuple containing two elements:
               1. A generator object that iterates over the connected components of the graph.
               2. A graph object representing the input list as a graph.
    """
    graph = nx.Graph()
    graph.add_nodes_from(list(range(0, number_of_atoms)))
    graph.add_edges_from(atom_pairs)
    return (nx.connected_components(graph), graph)


def create_atompair(matrix_contact) : 
    """Create a list of atom pairs based on the contact matrix.

    This function takes a contact matrix as input and returns a list of atom pairs
    where their value (their index in the distance matrix) is below a certain threshold.

    Parameters
    ----------
        matrix_contact : numpy.ndarray
            A 2D numpy array representing the contact matrix.

    Returns
    -------
        numpy.ndarray
            A 2D numpy array containing pairs of atoms where the value in the contact matrix
        is below the threshold. Each row represents an atom pair.
    """
    atom_pair = np.argwhere(matrix_contact)  # list of pair where their value (their indice in the distance matrix) is below the threshold
    atom_pair = atom_pair[atom_pair[:, 0] != atom_pair[:, 1]]
    return atom_pair


def relabel_node(graph, G, atom_names, resnames):
    """Iterate each connected subgraph, and relabel their node by [atom_name]_[index].

    Args:
        graph (generator): An object that lists all connected graphs
        G (networkx.classes.graph.Graph): The graph for this file
        atom_name (numpy.ndarray): A list of each atoms for this file
        resnames ():  

    Returns:
        list: list where each node is relabel
    """
    modified_subgraphs = []   # List to store modified subgraphs
    
    # For each connected component, create a subgraph and add it to the list
    for component in graph:
        # Create a subgraph from the component
        subgraph = G.subgraph(component)

        # Create a replacement dictionary for node names
        mapping_atom_name = {node: [atom_names[node], resnames[node]] for node in subgraph.nodes()}
        
        # Create a new subgraph with modified node names
        for node, (atom_name, resname) in mapping_atom_name.items():
            nx.set_node_attributes(subgraph, {node:node}, "label")
            nx.set_node_attributes(subgraph, {node:atom_name}, name="atom_name")
            nx.set_node_attributes(subgraph, {node:resname}, name="resnames")
        modified_subgraphs.append(subgraph)  # Add the modified subgraph to the list of modified subgraphs
    return modified_subgraphs


def get_graph_fingerprint(g):
    nodes = g.number_of_nodes()
    edges = g.number_of_edges()
    atom_names = " ".join(sorted(nx.get_node_attributes(g, "atom_name").values()))
    resnames = " ".join(sorted(set((nx.get_node_attributes(g, "resnames").values()))))
    return(nodes, edges, atom_names, resnames)


def print_groupby(object_groupby): 
    for f, g in object_groupby:
        print(f"graph_fingerprint {f[0]} | {f[1]} | {f[2]} | {f[3]}")
        for i in list(g): 
            print("\t", i.nodes()(data=True))
    print("\n")


def count_molecule (list_graph): 
    """_summary_
    https://stackoverflow.com/questions/46999771/comparing-a-large-number-of-graphs-for-isomorphism


    Args:
        list_graph (_type_): _description_

    Returns:
        _type_: _description_
    """
    dict_count = {}

    sorted_graphs = sorted(list_graph, key=get_graph_fingerprint)
    # print_groupby(groupby(sorted_graphs, key=get_graph_fingerprint))

    for fingerprint, graph in groupby(sorted_graphs, key=get_graph_fingerprint): 
        # fingerprint : (nb_node, nb_edge, atom_name)
        # graph : objet itertools that group all graph with the same fingerprint
        similar_graphs = list(graph)  # A list that contain all graph with the same fingerprint
        nb_graph = len(similar_graphs)  # Number of graph for this fingerprint

        if nb_graph == 1:  # For this fingerprint, there is only one graph
            dict_count[similar_graphs[0]] = nb_graph
        else: 
            if fingerprint[0] == 1:  # For this fingerprint, all the graph only have one node 
                dict_count[similar_graphs[0]] = nb_graph
            else:
                dict_count[similar_graphs[0]] = nb_graph
    return dict_count


def print_count(count, option=""):
    """Print the count of molecules.

    Args:
        count (dict): A dictionary containing information about the molecules.
        option (str): either we want to print the composition of each molecule (with the option "detail") or not
    """
    logger.info("Here is the content of this GRO file :")
    size = 0
    for i, (key, value) in enumerate(count.items()):
        if option=="detail":
            atom_names = " ".join(nx.get_node_attributes(key, "atom_name").values())
            logger.success(f"\nMolecule {i+1}: \n\t {key.number_of_nodes():,} atoms \n\t Quantity: {value:,} \n\t Composition: {atom_names}")
        else : 
            logger.success(f"\nMolecule {i+1}: \n\t {key.number_of_nodes():,} atoms \n\t Quantity: {value:,}")
        size += value
    logger.success(f"It containt {size:,} molecules")


def print_graph(dict_graph_count):
    for graph in dict_graph_count.keys():
        plt.figure()
        nx.draw(graph, node_color="green", with_labels=True, labels=nx.get_node_attributes(graph, "atom_name"))
        plt.show()


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
    logger.success(f"Threshold: {threshold} Angstrom")

    file_gro = mda.Universe(filepath_gro)  # load the .gro file
    logger.success(f"How many atoms there is in this file (at the beginning): {len(file_gro.atoms):,}")

    file_gro = delete_hydrogen_grofile(filepath_gro)
    logger.success(f"How many atoms there is in this file (without all hydrogen): {len(file_gro.atoms):,}")

    # https://docs.mdanalysis.org/1.1.0/documentation_pages/analysis/distances.html#MDAnalysis.analysis.distances.contact_matrix
    logger.info("Create contact_matrix ...")
    matrix_contact = contact_matrix(file_gro.atoms.positions, cutoff=threshold, returntype='sparse')

    # https://numpy.org/doc/stable/reference/generated/numpy.argwhere.html
    # https://www.includehelp.com/python/how-to-get-indices-of-elements-that-are-greater-than-a-threshold-in-2d-numpy-array.aspx
    logger.info("Create list of pair ...")
    atom_pair = create_atompair(matrix_contact)

    # Turn each tuple into a graph, if one node's label is already existant
    # It attached the other label to it
    # So it connected all the node who have common label
    logger.info("Convert atom pairs to graph ...")
    connexgraph_return, graph_return = convert_atompairs_to_graph(atom_pair, len(file_gro.atoms))

    logger.info("Begin relabel_node ...")
    list_graph = relabel_node(connexgraph_return, graph_return, file_gro.atoms.names, file_gro.resnames)

    logger.info("Counting molecules version1...")
    dict_count = count_molecule(list_graph)

    logger.info("Print dictionnary count and graph...")
    print_count(dict_count)
    # print_graph(dict_count)
    
    print("[main] ------------")

