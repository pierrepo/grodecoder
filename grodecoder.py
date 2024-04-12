"""Extract each molecule of a GRO file and print their occurence.

Usage:
    import grodecoder as gd
"""

__authors__ = "Karine DUONG"
__contact__ = "karine.duong@etu.u-paris.fr"
__copyright__ = "IBPC"
__date__ = "2024-03-18"


from collections import Counter
from itertools import groupby
from pathlib import Path
import re

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


def get_distance_matrix_between_atom(file_gro):
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


def add_attributes_to_nodes(graph, mol_system):
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
            The input graph.

    Returns
    -------
        list
            A list of subgraphs, each subgraph representing a connected component of the input graph.
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
    res_names = " ".join(sorted(set((nx.get_node_attributes(graph, "residue_name").values()))))

    graph_degrees = Counter(dict(graph.degree).values())
    degree_dist = " ".join([f"{key}:{value}" for key, value in sorted(graph_degrees.items())])
    return (nodes, edges, atom_names, res_names, degree_dist)


def print_graph_fingerprint(graph):
    """Print a graph fingerprint.

    Parameters
    ----------
        graph : networkx.classes.graph.Graph
            A NetworkX graph object.
    """
    logger.debug("print groupby ... ")
    fingerprint = get_graph_fingerprint(graph)
    logger.debug("Graph fingerprint-----------------")
    logger.debug(f"- Number of nodes: {fingerprint[0]}")
    logger.debug(f"- Number of edges: {fingerprint[1]}")
    logger.debug(f"- Sorted atom names (first 50 char.): {fingerprint[2][:50]}")
    logger.debug(f"- Sorted set of residue names: {fingerprint[3]}")
    logger.debug(f"- Node degrees dist: {fingerprint[4]}")


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


def print_graph_inventory(graph_dict):
    """Print graph inventory.

    Parameters
    ----------
        graph_dict: dict
            A dictionary with graphs as keys and counts (numbers of graphs) as values.
    """
    logger.info("Molecules inventory:")
    total_molecules_count = 0
    for graph_idx, (graph, count) in enumerate(graph_dict.items(), start=1):
        logger.info(f"Molecule {graph_idx:,} ----------------")
        logger.info(f"- number of atoms: {graph.number_of_nodes():,}")
        logger.info(f"- number of molecules: {count:,}")
        atom_names = list(nx.get_node_attributes(graph, "atom_name").values())
        atom_names_str = " ".join( atom_names[:20] )
        logger.debug(f"- 20 first atom names: {atom_names_str}")
        total_molecules_count += count
    logger.success(f"{total_molecules_count:,} molecules in total")


def print_graph(graph, filepath_name, option_color=False):
    """Print and save a graph as PNG

    Ressources
    ----------
    - Ressources to explore to increase node spacing with networkx-spring-layout 
    https://stackoverflow.com/questions/14283341/how-to-increase-node-spacing-for-networkx-spring-layout

    Parameters
    ----------
        graph : networkx.classes.graph.Graph
            A NetworkX graph object.
        filepath_name : str 
            Filepath to where I want to save the output graph
        option_color: str
            Either we want to display the nodes of the graph colored. By default, it's False.

    """
    plt.figure()

    if option_color:
        node_colors = []
        for node in graph.nodes:
            atom_name = re.sub(r'\d', '', graph.nodes[node]['atom_name'])
            if atom_name == 'C' or atom_name == 'CA' or atom_name == 'CB' or atom_name == 'CD' or atom_name == 'CE' or atom_name == 'CG' or atom_name == 'CZ':
                node_colors.append('black')
            elif atom_name == 'O' or atom_name == 'OE' or atom_name == 'OH' or atom_name == 'OD' or atom_name == 'OT' or atom_name == 'OG':
                node_colors.append('red')
            elif atom_name == 'N' or atom_name == 'ND' or atom_name == 'NE' or atom_name == 'NH' or atom_name == 'NZ':
                node_colors.append('blue')
            else:
                node_colors.append('green')  # Default color for other labels

        nx.draw(graph, node_color=node_colors, node_size = 75,
                with_labels=True, labels=nx.get_node_attributes(graph, "atom_name"),
                edge_color = "grey")
    else:
        nx.draw(graph, node_color="green", node_size = 75,
            with_labels=True, labels=nx.get_node_attributes(graph, "atom_name"),
            edge_color = "grey")    
    plt.savefig(filepath_name)


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


def read_structure_file_remove_hydrogens(file_path):
    """Read Groamacs .gro file and remove hydrogen atoms.

    Parameters
    ----------
        filepath: str
            The filepath of the structure file (.gro, .pdb)..

    Returns
    -------
        MDAnalysis.core.universe.Universe
            A MDAnalysis universe object containing only non-hydrogen atoms.
    """
    logger.info(f"Reading file: {file_path}")
    molecule = mda.Universe(file_path)
    logger.success(f"Found {len(molecule.atoms):,} atoms")
    # Print 10 first atoms for debugging.  
    print_first_atoms(molecule)
    logger.info("Removing H atoms...")
    molecule_without_h = molecule.select_atoms("not (name H*)")
    logger.success(f"{len(molecule_without_h.atoms):,} atoms remaining")
    # Print 10 first atoms for debugging.  
    print_first_atoms(molecule_without_h)
    without_h_file_path = Path(file_path).stem + "_without_H.gro"
    molecule_without_h.write(without_h_file_path)
    return molecule_without_h


def check_overlapping_residue_between_graphs(graph_list):
    """Check there is no overlapping residue between graphs.

    This function checks that there are no overlapping residue between different graphs/molecules.
    It extracts  the set of residue IDs for each graph/molecule and performs an intersection 
    operation between these sets. If any intersection is found, it indicates
    that there are overlapping residue between graph/molecules.

    Parameters
    ----------
        graph_list: list
            A list of NetworkX graph objects representing different molecules.
    """
    logger.info("Verifying residue overlapping...")
    res_id_set_all = set()
    res_id_common = []

    for graph in graph_list:
        res_id_set = set((nx.get_node_attributes(graph, "residue_id").values()))
        res_id_intersect = res_id_set_all.intersection(res_id_set)
        res_id_set_all.update( res_id_set )
        if res_id_intersect:
            for res_id in res_id_intersect:
                logger.critical(f"Residue id {res_id} is found in multiple graphs")
                res_id_common.append(res_id)

    if not res_id_common:
        logger.success("No overlapping residue found")
    else:
        for res_id in res_id_common:
            for graph_id, graph in enumerate(graph_list):
                res_id_set = set((nx.get_node_attributes(graph, "residue_id").values()))
                if res_id in res_id_set:
                    logger.error(f"Residue id {res_id} is found in graph id {graph_id}")
        raise Exception("Some residue id are found in multiple graphs")


def is_protein(graph):
    """Check if the molecule represented by the graph is a protein.

    This function checks whether the graph represents a protein molecule by 
    verifying the presence of C-alpha (CA) atoms and their corresponding amino acids.

    Parameters
    ----------
        graph: networkx.Graph
            The input graph representing the molecule.

    Returns
    -------
        bool
            True if the molecule is a protein, False otherwise.
    """
    logger.info("Checking if the molecule is a protein...")
    atom_names = ""
    if graph.number_of_nodes() > 1:
        atom_names = " ".join(sorted(nx.get_node_attributes(graph, "atom_name").values()))
    return ("CA" in atom_names)


def extract_protein_sequence(graph):
    """Extract the protein sequence from a graph.

    This function extracts the protein sequence from the molecule represented 
    by the input graph. It looks for the residue names corresponding to the C-alpha 
    (CA) atoms in each node and converts them to single-letter amino acid codes.

    The correspondence between 3-letter code and 1-letter code
    is taken from MDAnalysis:
    https://docs.mdanalysis.org/1.0.1/_modules/MDAnalysis/lib/util.html#convert_aa_code
    See also the list of known residue names in MDAnalysis:
    https://userguide.mdanalysis.org/stable/standard_selections.html#proteins

    Parameters
    ----------
        graph: networkx.Graph
            The input graph representing the molecule.

    Returns
    -------
        str
            Protein sequence.
    """
    logger.info("Extractigg protein sequence...")
    AMINO_ACID_DICT = mda.lib.util.inverse_aa_codes
    protein_sequence = []
    # graph.nodes.items() returns an iterable of tuples (node_id, node_attributes).
    # Examples:
    # (0, {'atom_name': 'N', 'residue_id': 1, 'residue_name': 'MET', 'label': 1})
    # (1, {'atom_name': 'CA', 'residue_id': 1, 'residue_name': 'MET', 'label': 5})
    # (2, {'atom_name': 'CB', 'residue_id': 1, 'residue_name': 'MET', 'label': 7})
    for _, node_attr in sorted(graph.nodes.items(), key=lambda x: x[0]):
        if node_attr["atom_name"] == "CA":
            protein_sequence.append( AMINO_ACID_DICT.get(node_attr["residue_name"], "?") )
    return "".join(protein_sequence)


def main(input_file_path, draw_graph_option=False, check_overlapping_residue=False):
    """Excute the main function for analyzing a .gro file.

    Parameters
    ----------
        input_file_path: str
            Filepath of the .gro file we want to analyzed
        draw_graph_option: boolean
            Draw the graph of each molecule and save it as a PNG file. Default: False.
        check_overlapping_residue: boolean
            Check of some residues are overlapping between graphs / molecules. Default: False.
    """
    threshold = max(BOND_LENGTH.values())
    logger.success(f"Bond threshold: {threshold} Angstrom")

    molecular_system = read_structure_file_remove_hydrogens(input_file_path)

    atom_pairs = get_atom_pairs(molecular_system, threshold)

    graph_return = convert_atom_pairs_to_graph(atom_pairs, len(molecular_system.atoms))

    graph_with_node_attributes = add_attributes_to_nodes(graph_return, molecular_system)
    graph_list = get_graph_components(graph_with_node_attributes)

    if check_overlapping_residue:
        check_overlapping_residue_between_graphs(graph_list)

    graph_count_dict = count_molecule(graph_list)

    # Print fingerprint for each graph/molecule
    for graph in graph_count_dict.keys():
        print_graph_fingerprint(graph)

    print_graph_inventory(graph_count_dict)

    if draw_graph_option:
        logger.info("Drawing graphs...")
        filename = Path(input_file_path).stem
        for index_graph, graph_count in enumerate(graph_count_dict.keys()):
            print_graph(graph_count, f"{filename}_{index_graph}.png")

    protein_sequence = {}
    for index_graph, graph in enumerate(graph_count_dict.keys()):
        if is_protein(graph):
            protein_sequence[index_graph] = extract_protein_sequence(graph)
    
    logger.info("Print the protein sequence...")
    for index_graph, seq in protein_sequence.items():
        print(f"{index_graph}: {seq}")


def is_a_structure_file(filepath):
    """Check if the given filepath points to an existing structure file (.gro, .pdb).

    Parameters
    ----------
        filepath : str
            Path of the file.

    Raises
    ------
        argparse.ArgumentTypeError
            If the given filepath is not an existing file,
            or if it does not have a '.gro' or '.pdb' extension.

    Returns
    -------
        str
            The validated path.
    """
    filename = Path(filepath)
    if not Path.is_file(filename):
        raise argparse.ArgumentTypeError(f"{filepath} does not exist")

    if filename.suffix not in [".gro", ".pdb"]:
        raise argparse.ArgumentTypeError(f"{filepath} is not a .gro or .pdb file.")
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
                                     description="Extract molecules from a structure file (.gro, .pdb).",
                                     usage="grodecoder.py [-h] --input structure_file [--drawgraph]")
    parser.add_argument("--input",
                        type=is_a_structure_file,
                        help="structure file path (.gro, .pdb)",
                        required=True)
    parser.add_argument("--drawgraph",
                        help="Draw graph of each molecule. Default: False.",
                        default=False,
                        action="store_true")
    parser.add_argument("--checkoverlapping",
                        help="Check if some residues are overlapping between residues. Default: False.",
                        default=False,
                        action="store_true")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_arg()
    main(args.input, draw_graph_option=args.drawgraph, check_overlapping_residue=args.checkoverlapping)