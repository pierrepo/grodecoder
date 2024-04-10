"""Extract each molecule of a GRO file and print their occurence.

Usage:
    import grodecoder as gd
"""

__authors__ = "Karine DUONG"
__contact__ = "karine.duong@etu.u-paris.fr"
__copyright__ = "IBPC"
__date__ = "2024-03-18"


from itertools import groupby
from pathlib import Path

import argparse
from loguru import logger
import MDAnalysis as mda
from MDAnalysis.analysis.distances import contact_matrix
import matplotlib.pyplot as plt
import numpy as np
import networkx as nx
from networkx.algorithms.components.connected import connected_components
from scipy.sparse import triu
from scipy.spatial.distance import cdist


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

    Warning
    -------
    MDAnalysis convert the coordonate in Angstrom
    despite in the gro file it's in nm
    https://userguide.mdanalysis.org/stable/units.html#table-baseunits
    https://manual.gromacs.org/current/reference-manual/file-formats.html#gro

    Parameters
    ----------
        file_gro: MDAnalysis.core.universe.Universe
            An object representing the molecular structure loaded from a GRO file.

    Returns
    -------
        numpy.ndarray
            matrix of interatomic distances between all atoms
    """
    position = file_gro.atoms.positions
    taille = len(file_gro.atoms)
    tmp = np.full((taille, taille), 100.0)
    for index_i, pos_i in enumerate(position):
        for index_j, pos_j in enumerate(position[index_i+1:], start=index_j+1):
            tmp[index_i][index_j] = cdist([pos_i], [pos_j])
    return tmp


def convert_atom_pairs_to_graph(atom_pairs, total_number_of_atoms):
    """Convert a list of pairs to a graph and its connected components.

    Reference
    ---------
    - https://stackoverflow.com/questions/4842613/merge-lists-that-share-common-elements

    Parameters
    ----------
        atom_pairs: list
            A list of pairs representing edges in the graph.
        number_of_atoms: int
            The total number of atoms in the graph.

    Returns
    -------
        tuple
            A tuple containing two elements:
               1. A generator object that iterates over the connected components of the graph.
               2. A graph object representing the molecular system.
    """
    logger.info("Converting atom pairs to graph...")
    graph = nx.Graph()
    # Add all atoms as single nodes.
    graph.add_nodes_from(list(range(0, total_number_of_atoms)))
    # Add atom pairs as edges.
    graph.add_edges_from(atom_pairs)
    return graph


def get_atom_pairs(molecular_system, threshold):
    """Create a list of atom pairs based on the contact matrix (which based on the input system and threshold).

    This function calculates the pairs of atoms in the molecular system based on their distances.
    It uses a specified threshold to determine which pairs of atoms are considered to be in contact.

    Reference
    ---------
    - https://docs.mdanalysis.org/1.1.0/documentation_pages/analysis/distances.html#MDAnalysis.analysis.distances.contact_matrix
    - https://numpy.org/doc/stable/reference/generated/numpy.argwhere.html

    Parameters
    ----------
        molecular_system: MDAnalysis.core.groups.AtomGroup
            The molecular system object containing information about atoms.
        threshold: float
            The distance threshold used to determine atom contacts.

    Returns
    -------
        numpy.ndarray
            An array containing pairs of atom indices representing atom contacts.
    """
    logger.info("Creating contact_matrix...")
    matrix = contact_matrix(molecular_system.atoms.positions, cutoff=threshold, returntype="sparse")
    # Output matrix is sparse matrix of type: scipy.sparse.lil_matrix
    # Keep only the upper triangular part of the sparse matrix (without the diagonal)
    # https://docs.scipy.org/doc/scipy-1.13.0/reference/generated/scipy.sparse.triu.html
    matrix = triu(matrix, k=1)

    logger.info("Creating atom pairs list...")
    atom_pairs = np.argwhere(matrix)
    logger.success(f"Found {len(atom_pairs):,} atom pairs")
    return atom_pairs


def relabel_node(graph, mol_system):
    """Iterate each connected subgraph, and relabel their node by [atom_name]_[index].

    Parameters
    ----------
        old_graph: networkx.classes.graph.Graph
            The graph for this file
        mol_system: MDAnalysis.core.groups.AtomGroup
            The molecular system object containing information about atoms.
            
    Returns
    -------
        list
            list where each node is relabel
    """
    logger.info("Relabeling nodes in graph...") 
    logger.success(f"Old graph: {graph.number_of_nodes():,} nodes")

    # Create a replacement dictionary for the new attributes for each node
    atoms_to_matrix_ids = {}
    for node_id in graph.nodes():
        atoms = {}
        atoms["atom_id"] = mol_system.atoms.ids[node_id]
        atoms["atom_name"] = mol_system.atoms.names[node_id]
        atoms["res_id"] = mol_system.resids[node_id]
        atoms["res_name"] = mol_system.resnames[node_id]

        atoms_to_matrix_ids[node_id] = atoms

    for node_id, atom in atoms_to_matrix_ids.items():
        # Change properties first.
        nx.set_node_attributes(graph, {node_id: atom["atom_name"]}, name="atom_name")
        nx.set_node_attributes(graph, {node_id: atom["res_id"]}, name="residue_id")
        nx.set_node_attributes(graph, {node_id: atom["res_name"]}, name="residue_name")
        # Then node id.
        nx.set_node_attributes(graph, {node_id: atom["atom_id"]}, "label")
    logger.success(f"New graph: {graph.number_of_nodes():,} nodes")
    return graph


def get_graph_components(graph):
    """Extract the connected components of a graph.

    Parameters
    ----------
        graph: networkx.Graph
            The input graph from which to extract connected components.

    Returns
    -------
        list
            A list containing subgraphs, each representing a connected component of the input graph.
    """
    logger.info("Extracting graph components...")
    graph_components = connected_components(graph)
    graph_list = [graph.subgraph(subgraph) for subgraph in graph_components]

    logger.success(f"Found {len(graph_list):,} subgraphs")
    return graph_list


def get_graph_fingerprint(graph):
    """Generate a fingerprint for a given graph.

    This function calculates a fingerprint for a given graph based on its properties, including
    the number of nodes, the number of edges, and sorted concatenations of atom names, residue 
    names and their degree.
    
    Reference
    ---------
    - https://stackoverflow.com/questions/46999771/comparing-a-large-number-of-graphs-for-isomorphism

    Parameters
    ----------
        graph: networkx.classes.graph.Graph
            The graph for which the fingerprint is to be generated.

    Returns
    -------
        tuple:
            A tuple containing the fingerprint of the graph, which includes the following elements:
               - Number of nodes in the graph.
               - Number of edges in the graph.
               - Concatenation of sorted atom names of all nodes in the graph.
               - Concatenation of sorted residue names of all nodes in the graph.
               - Concatenation of sorted degree of all nodes in the graph.
    """
    nodes = graph.number_of_nodes()
    edges = graph.number_of_edges()
    atom_names = " ".join(sorted(nx.get_node_attributes(graph, "atom_name").values()))
    resnames = " ".join(sorted(set((nx.get_node_attributes(graph, "resnames").values()))))

    dict_degree = {index: value for index, (key, value) in enumerate(sorted(graph.degree))}
    degree = " ".join([f"{key}:{value}" for key, value in dict_degree.items()])
    return (nodes, edges, atom_names, resnames, degree)


def print_groupby(object_groupby):
    """Print the grouped graphs along with their fingerprints.

    This function iterates over the groups created by the input groupby object and prints
    the fingerprint along with the nodes and their attributes for each graph in the group.

    Parameters
    ----------
        object_groupby: itertools.groupby
            An object generated by the groupby function.
    """
    logger.info("print groupby ... ")
    for fingerprint, graph in object_groupby:
        print(f"graph_fingerprint {fingerprint[0]} | {fingerprint[1]} | {fingerprint[2]} | {fingerprint[3]} | {fingerprint[4]}")
        for subgraph in list(graph):
            print("\t", subgraph.nodes()(data=True))
    print("\n")


def count_molecule(graph_list):
    """Count the occurrence of molecules in a list of graphs based on their fingerprints.

    This function takes a list of graphs and counts the occurrence of each unique molecule
    based on their fingerprints, which are calculated using the get_graph_fingerprint function.
    
    Reference 
    ---------
    - https://stackoverflow.com/questions/46999771/comparing-a-large-number-of-graphs-for-isomorphism

    Parameters
    ----------
        graph_list: list
            A list of graph objects.

    Returns
    -------
        dict
            A dictionary where keys are unique molecules (graph objects) and values are\
                their respective occurrence counts in the input list of graphs.
    """
    logger.info("Counting molecules...")
    dict_count = {}
    sorted_graphs = sorted(graph_list, key=get_graph_fingerprint)
    # print_groupby(groupby(sorted_graphs, key=get_graph_fingerprint))

    for fingerprint, graph in groupby(sorted_graphs, key=get_graph_fingerprint):
        # fingerprint : (nb_node, nb_edge, atom_name, resname, degree)
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


def print_count(count, option):
    """Print the count of molecules.

    Parameters
    ----------
        count: dict
            A dictionary containing information about the molecules.
        option: str
            Either we want to print the composition of each molecule (with the option True) or not
    """
    logger.info("File content:")
    total_molecules_count = 0
    for mol_idx, (mol_graph, mol_count) in enumerate(count.items(), start=1):
        logger.success(f"Molecule {mol_idx:,} ----------------")
        logger.success(f"- number of atoms: {mol_graph.number_of_nodes():,}")
        logger.success(f"- number of molecules: {mol_count:,}")
        if option:
            atom_names = " ".join(nx.get_node_attributes(mol_graph, "atom_name").values())
            logger.success(f"- atom names: {atom_names}")
        total_molecules_count += mol_count
    logger.success(f"{total_molecules_count:,} molecules in total")


def print_graph(dict_graph_count):
    """Display each graph in a dictionary along with its count.

    This function takes a dictionary where the keys are graph objects and the values
    are their respective counts. It then displays each graph using Matplotlib's pyplot
    library, with the graph nodes colored green and labeled with their atom names.

    Parameters
    ----------
        dict_graph_count: dict
            A dictionary where the keys are graph objects and the values are their respective counts.
    """
    for graph in dict_graph_count.keys():
        plt.figure()
        nx.draw(graph, node_color="green", with_labels=True, labels=nx.get_node_attributes(graph, "atom_name"))
        plt.show()


def print_first_atoms(mda_universe, number_of_atoms=10):
    """Print the first atoms in the MDAnalysis Universe object.

    For debugging purpose only.

    Parameters
    ----------
        mda_universe: MDAnalysis.core.universe.Universe
            An object representing the molecular structure.
        number_of_atoms: int
            The number of atoms to be printed. Default is 10.   
    """
    for atom in mda_universe.atoms[:number_of_atoms]:
        string_coords = " ".join([f"{coord:.2f}" for coord in atom.position])
        logger.debug(f"Res. id: {atom.residue.resid} | res. name: {atom.residue.resname} | "
                     f"at. name: {atom.name} | at. id: {atom.id} | "
                     f"at. coord.: {string_coords}"
        )


def read_gro_files_remove_hydrogens(gro_file_path):
    """Read Groamacs .gro file and remove hydrogen atoms.

    Parameters
    ----------
        filepath: str
            The filepath of the .gro file containing hydrogen atoms.

    Returns
    -------
        molecule_without_H
            A modified MDA Universe object containing only non-hydrogen atoms.
    """
    logger.info(f"Reading file: {gro_file_path}")
    molecule = mda.Universe(gro_file_path)
    logger.success(f"Found {len(molecule.atoms):,} atoms")
    # Print 10 first atoms for debugging.  
    # print_first_atoms(molecule)
    logger.info("Removing H atoms...")
    molecule_without_h = molecule.select_atoms("not (name H*)")
    logger.success(f"{len(molecule_without_h.atoms):,} atoms remaining")
    # Print 10 first atoms for debugging.  
    # print_first_atoms(molecule_without_h)
    without_h_file_path = "./results/" + Path(gro_file_path).stem + "_withoutH.gro"
    molecule_without_h.write(without_h_file_path)
    return molecule_without_h


def control_quality (graph_list):
    """Check the quality of the molecular graphs.

    This function checks the quality of the molecular graphs to ensure that 
    there are no overlapping atoms between different molecules. It extracts 
    the set of residue IDs for each molecule and performs an intersection 
    operation between these sets. If any intersection is found, it indicates
    that there are overlapping atoms between molecules.

    Parameters
    ----------
        graph_list: list
            A list of NetworkX graph objects representing different molecules.
    """
    logger.info("Quality control ...")
    list_set = []
    for component in graph_list:
        atom_id = set((nx.get_node_attributes(component, "res_id").values()))
        list_set.append(atom_id)

    for index_i in range(len(list_set)-1):
        for index_j in range(index_i+1, len(list_set)):
            intersect = list_set[index_i].intersection(list_set[index_j])
            if intersect != set():
                logger.warning(f"Intersection between molecules {index_i} and {index_j} with the atom {intersect}")
    logger.success("No intersection between atom of different molecule")


def main(filepath_gro, print_molecule_option, print_graph_option):
    """Excute the main function for analyzing a .gro file.

    Parameters
    ----------
        filepath_gro: str
            Filepath of the .gro file we want to analyzed
        print_molecule_option: boolean
            Either we want to print the composition of each molecule (with the option True) or not
        print_graph_option: boolean
            Either we want to print the graph of each molecule (with the option True) or not
    """
    threshold = max(BOND_LENGTH.values())
    logger.success(f"Threshold: {threshold} Angstrom")

    molecular_system = read_gro_files_remove_hydrogens(filepath_gro)

    atom_pair = get_atom_pairs(molecular_system, threshold)

    graph_return = convert_atom_pairs_to_graph(atom_pair, len(molecular_system.atoms))

    graph_relabel = relabel_node(graph_return, molecular_system)
    graph_connex_list = get_graph_components(graph_relabel)

    dict_count = count_molecule(graph_connex_list)

    control_quality(graph_connex_list)

    logger.info("Print dictionnary count and graph...")
    print_count(dict_count, print_molecule_option)
    if print_graph_option:
        print_graph(dict_count)

    print("---------------")


def is_an_existing_gro_file(filepath):
    """Check if the given filepath points to an existing GRO file.

    Parameters
    ----------
        filepath : str
            The path to be checked.

    Raises
    ------
        argparse.ArgumentTypeError
            If the given filepath is not a file or does not exist, or if it does not have '.gro' extension

    Returns
    -------
        str
            The validated path.
    """
    source = Path(filepath)
    if not Path.is_file(source):
        raise argparse.ArgumentTypeError(f"{filepath} not exist")

    if Path(filepath).suffix != ".gro":
        raise argparse.ArgumentTypeError(f"{filepath} is not a GRO file.")
    return filepath


def parse_arg():
    """Parse command-line arguments.

    This function uses the argparse module to parse the command-line arguments
    provided by the user. It sets up the argument parser with information about
    the program, its usage, and the available options.

    Ressources
    ----------
    - https://docs.python.org/3/library/argparse.html

    Return
    ------
        argparse.Namespace
            An object containing the parsed arguments as attributes.

    """
    parser = argparse.ArgumentParser(prog="grodecoder",
                                     description="Programm to extract each molecule of a GRO file and print their occurence.",
                                     usage="grodecoder.py [-h] -g GRO [-pm PRINTMOLECULE] [-pg PRINTGRAPH]")

    parser.add_argument("-g", "--gro",
                        type=is_an_existing_gro_file,
                        help="a GRO filepath in input to this programm",
                        required=True)

    parser.add_argument('-pm', "--printmolecule",
                        help="Either we want to print the composition of each molecule (with the option True) or not. By default it's False.",
                        default=False)

    parser.add_argument('-pg', "--printgraph",
                        help="Either we want to print the graph of each molecule (with the option True) or not. By default it's False.",
                        default=False)
    return parser.parse_args()


if __name__=="__main__":
    args = parse_arg()
    main(args.gro, args.printmolecule, args.printgraph)
