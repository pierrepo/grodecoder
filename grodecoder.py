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
from collections import Counter


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
        mapping_atom_name = {node: atom_name[node] for node in subgraph.nodes()}
        # Create a new subgraph with modified node names
        for node, label in mapping_atom_name.items():
            nx.set_node_attributes(subgraph, {node:node}, "label")
            nx.set_node_attributes(subgraph, {node:label}, name="atom_name")
        modified_subgraphs.append(subgraph)  # Add the modified subgraph to the list of modified subgraphs
    return modified_subgraphs


def compare_nodes(G1, G2):
    return G1["atom_name"] == G2["atom_name"]


def sort_by_size (list_graph) : 
    dict_size = {}  # Dictionary to store each subgraph based on their size
    for i in range(len(list_graph)): # Iterate on each sublist, and sort it based on their size
        size = list_graph[i].number_of_nodes()
        if size in dict_size:
            dict_size[size] += [i]
        else:
            dict_size[size] = [i]

    dict_size = dict(sorted(dict_size.items()))
    for key, value in dict_size.items():
        logger.info(f"{len(value)} graph with {key} nodes")

    return dict_size


def count_molecule(list_graph, dict_size):
    """Print information about the modified subgraphs.

    Args:
        list_graph (list): A list of modified subgraphs.
        dict_size (dict): A dictionary with in key the size of a graph, and in value all the index of graph (from list_graph) that have this size

    Returns:
        dict: A dictionary containing information about the subgraphs.
              Keys are tuples representing the node labels of each subgraph,
              and values are the number of occurrences of each subgraph.
    """
    dict_count = {}  # Dictionary to store each subgraph based on their occurence

    for taille, subgraph_list in dict_size.items():  # Iterate on the dictionary
        k = 0
        while k < len(subgraph_list):
            if len(subgraph_list) == 1:  # If the number of subgraph for this size is 1
                name_molecule_list1 = tuple(name for name in nx.get_node_attributes(subgraph_list[k], "atom_name").values())
                dict_count[name_molecule_list1] = 1  # Add it to the dictionary  
            else : 
                if taille==1 :  # See if the subgraph only have one node
                    for name in nx.get_node_attributes(subgraph_list[k], "atom_name").values():
                        dict_count[name] = dict_count.get(name, 0) + 1 
                else : 
                    i = k 
                    while i < len(subgraph_list)-1 : 
                        if nx.is_isomorphic(subgraph_list[i], subgraph_list[i+1], node_match=compare_nodes): 
                            name_molecule_list1 = tuple(sorted(name for name in nx.get_node_attributes(subgraph_list[i], "atom_name").values()))
                            dict_count[name_molecule_list1] = dict_count.get(name_molecule_list1, 1) + 1 
                        
                        else : 
                            name_molecule_list1 = tuple(sorted(name for name in nx.get_node_attributes(subgraph_list[i], "atom_name").values()))
                            name_molecule_list2 = tuple(sorted(name for name in nx.get_node_attributes(subgraph_list[i+1], "atom_name").values()))
                             
                            dict_count[name_molecule_list1] = dict_count.get(name_molecule_list1, 0) + 1
                            dict_count[name_molecule_list2] = dict_count.get(name_molecule_list2, 0) + 1
                        i += 1 
                        k = i+1
            k += 1
    return dict_count


def count_molecule2(list_graph, dict_size):
    dict_count = {}  # "atom_name" : occurence 
    dict_graph = {}  # "atom_name" : graph

    for nb_node, subgraph_index_list in dict_size.items():
        if nb_node == 1:  # je n'ai qu'un node
            name_node = {index : list(nx.get_node_attributes(list_graph[index], "atom_name").values())[0] 
                         for index in subgraph_index_list
                         }
            occurence_name_node = Counter(name_node.values())
            dict_count = {**dict_count, **occurence_name_node}
            dict_graph = {value: key for key, value in name_node.items()}

        elif len(subgraph_index_list) == 1: 
            name_molecule = tuple(sorted(nx.get_node_attributes(list_graph[subgraph_index_list[0]], "atom_name").values()))
            dict_count[name_molecule] = 1
            dict_graph[name_molecule] = subgraph_index_list[0]

        else: 
            i = 0 
            while i < len(subgraph_index_list)-1: 
                index1 = subgraph_index_list[i] 
                index2 = subgraph_index_list[i+1]
                if nx.is_isomorphic(list_graph[index1], list_graph[index2], node_match=compare_nodes): 
                    name_molecule_list1 = tuple(sorted(nx.get_node_attributes(list_graph[index1], "atom_name").values()))
                    dict_count[name_molecule_list1] = dict_count.get(name_molecule_list1, 1) + 1
                    dict_graph[name_molecule_list1] = dict_graph.get(name_molecule_list1, index1)
                else:
                    name_molecule_list1 = tuple(sorted(nx.get_node_attributes(list_graph[index1], "atom_name").values()))
                    name_molecule_list2 = tuple(sorted(nx.get_node_attributes(list_graph[index2], "atom_name").values()))
                    dict_count[name_molecule_list1] = dict_count.get(name_molecule_list1, 0) + 1
                    dict_count[name_molecule_list2] = dict_count.get(name_molecule_list2, 0) + 1
                    dict_graph[name_molecule_list1] = dict_graph.get(name_molecule_list1, index1)
                    dict_graph[name_molecule_list2] = dict_graph.get(name_molecule_list2, index2)
                i +=1
    return (dict_count,dict_graph)


def get_graph_fingerprint(g):
    nodes = g.number_of_nodes()
    edges = g.number_of_edges()
    #atom_names = " ".join(sorted(set(nx.get_node_attributes(g, "atom_name").values())))  # pourquoi mettre des sets ? 
    atom_names = " ".join(sorted(nx.get_node_attributes(g, "atom_name").values()))
    return(nodes, edges, atom_names)


def print_groupby(object_groupby): 
    for f, g in object_groupby:
        print(f"graph_fingerprint {f[0]} | {f[1]} | {f[2]}")
        for i in list(g): 
            print("\t", i.nodes()(data=True))
    print("\n")


# https://stackoverflow.com/questions/46999771/comparing-a-large-number-of-graphs-for-isomorphism
from itertools import groupby 
def count_molecule3 (list_graph): 
    dict_count = {}
    dict_graph = {}
    
    sorted_graphs = sorted(list_graph, key=get_graph_fingerprint)
    # print_groupby(groupby(sorted_graphs, key=get_graph_fingerprint))

    for f, g in groupby(sorted_graphs, key=get_graph_fingerprint): 
        # f : (nb_node, nb_edge, atom_name)
        # g : objet itertools qui regroupe les graphes avec les memes caractéristiques f
        similar_graphs = list(g)  # [tous les graph qui ont les memes caractéristiques f]
        nb_graph = len(similar_graphs) #nb de graph pour chaque caractéristique
        atom_name = f[2]

        if nb_graph > 1:
            if f[0]==1:  # je n'ai qu'un node pour ces carac 
                for i in range(nb_graph):
                    dict_count[atom_name] = dict_count.get(atom_name, 0) + 1
                    dict_graph[atom_name] = dict_graph.get(atom_name, similar_graphs[i])
            else:
                for i in range(nb_graph):
                    for j in range(i + 1, nb_graph):
                        g1, g2 = similar_graphs[i], similar_graphs[j]
                        if g1 != g2 and nx.is_isomorphic(g1, g2):
                            dict_count[atom_name] = dict_count.get(atom_name, 1) + 1 
                            dict_graph[atom_name] = dict_graph.get(atom_name, g1)
                        else : 
                            dict_count[atom_name] = dict_count.get(atom_name, 0) + 1 
                            dict_count[atom_name] = dict_count.get(atom_name, 0) + 1 
                            dict_graph[atom_name] = dict_graph.get(atom_name, g1)
                            dict_graph[atom_name] = dict_graph.get(atom_name, g2)
        else:
            dict_count[atom_name] = dict_count.get(atom_name, 0) + 1 
            dict_graph[atom_name] = dict_graph.get(atom_name, similar_graphs[0])
    return (dict_count, dict_graph)


def reorganize_dictionnary(list_graph, dict_count, dict_graph) : 
    dict_return = {}
    for key, value in dict_count.items():
        value_graph = dict_graph[key]
        dict_return[list_graph[value_graph]] = value
    return dict_return


def reorganize_dictionnary3 (list_graph, dict_count, dict_graph): 
    dict_return = {}
    for key, value in dict_count.items():
        graph = dict_graph[key]
        dict_return[graph] = value
    return dict_return


def print_dict_graph_count (dict_graph_count, option=""):
    size = 0
    for i, (key, value) in enumerate(dict_graph_count.items()): 
        print(f"\nMolecule {i+1}: \n\t {len(key):,} atoms \n\t Quantity: {value:,}")
        if option=="detail":
            print(f"\t Composition: {key}")
        size += value
    print(f"It containt {size:,} molecules")


def print_graph(dict_graph_count):
    for graph in dict_graph_count.keys() : 
        nx.draw(graph, node_color="green", with_labels=True, labels=nx.get_node_attributes(graph, "atom_name"))





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
    logger.info("Create contact_matrix ...")
    mat_contact = contact_matrix(file_gro.atoms.positions, cutoff=threshold, returntype='sparse')

    # https://numpy.org/doc/stable/reference/generated/numpy.argwhere.html
    # https://www.includehelp.com/python/how-to-get-indices-of-elements-that-are-greater-than-a-threshold-in-2d-numpy-array.aspx
    logger.info("Create list of pair ...")
    pair_matrix = np.argwhere(mat_contact)  # list of pair where their value (their indice in the distance matrix) is below the threshold

    # Turn each tuple into a graph, if one node's label is already existant
    # It attached the other label to it
    # So it connected all the node who have common label
    logger.info("Create pairmatrix_to_graph ...")
    connexgraph_return, graph_return = pairmatrix_to_graph(pair_matrix)

    logger.info("Begin relabel_node ...")
    list_graph = relabel_node(connexgraph_return, graph_return, file_gro.atoms.names)

    logger.info("Begin sort_by_size()...")
    dict_sizegraph = sort_by_size (list_graph)

    logger.info("Counting molecules version1...")
    dict_count, dict_graph = count_molecule2(list_graph, dict_sizegraph)
    dict_graph_count = reorganize_dictionnary(list_graph, dict_count, dict_graph)

    # logger.info("Counting molecules version2...")
    # dict_count, dict_graph = count_molecule3(list_graph)
    # dict_graph_count = reorganize_dictionnary3(list_graph, dict_count, dict_graph)

    logger.info("Print dictionnary count and graph...")
    print_dict_graph_count(dict_graph_count)
    #print_graph(dict_graph_count) #--> pas final
    #print_count(count)
    
    print("[main] ------------")

