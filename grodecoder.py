"""Extract each molecule of a GRO file and print their occurence.

Usage:
    python grodecoder.py --input [input file] [other option]
"""

__authors__ = ("Karine DUONG", "Pierre POULAIN")
__contact__ = "pierre.poulain@u-paris.fr"

from collections import Counter
import datetime
import hashlib
import itertools
from itertools import groupby
import json
import os
import pandas as pd
from pathlib import Path
import subprocess
import time

import argparse
from loguru import logger
import MDAnalysis as mda
from MDAnalysis.analysis.distances import contact_matrix, self_distance_array
import numpy as np
import networkx as nx
from networkx.algorithms.components.connected import connected_components
from scipy.sparse import triu
from scipy.spatial.distance import cdist

import mol_def
import search_into_PDB


current_dir = os.path.dirname(os.path.abspath(__file__))

filepath_CSML = os.path.join(current_dir, "data/databases/lipid_CHARMM_GUI_CSML.csv")
CSML_CHARMM_GUI = pd.read_csv(filepath_CSML, sep=";")

filepath_MAD = os.path.join(current_dir, "data/databases/lipid_MAD.csv")
MAD_DB = pd.read_csv(filepath_MAD, sep=";")


def get_distance_matrix_between_atom(
    file_gro: mda.core.universe.Universe,
) -> np.ndarray:
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
    size = len(file_gro.atoms)
    tmp = np.full((size, size), 100.0)
    for index_i, pos_i in enumerate(position):
        for index_j, pos_j in enumerate(position[index_i + 1 :], start=index_i + 1):
            tmp[index_i][index_j] = cdist([pos_i], [pos_j])
    return tmp


def get_atom_pairs(
    molecular_system: mda.core.universe.Universe, threshold: float
) -> np.ndarray:
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
    matrix = contact_matrix(
        molecular_system.atoms.positions, cutoff=threshold, returntype="sparse"
    )
    # Output matrix is sparse matrix of type: scipy.sparse.lil_matrix
    # Keep only the upper triangular part of the sparse matrix (without the diagonal)
    # https://docs.scipy.org/doc/scipy-1.13.0/reference/generated/scipy.sparse.triu.html
    matrix = triu(matrix, k=1)

    logger.info("Creating atom pairs list...")
    atom_pairs = np.argwhere(matrix)
    logger.success(f"Found {len(atom_pairs):,} atom pairs")
    return atom_pairs


def get_atom_pairs_from_threshold(
    mol: mda.core.universe.Universe, threshold: float
) -> np.ndarray[np.ndarray]:
    """Get atom pairs within a specified distance threshold.

    Parameters
    ----------
        mol : MDAnalysis.core.universe.Universe
            The MDAnalysis universe object representing the molecular system.
        threshold : float
            The distance threshold for determining atom contacts.

    Returns
    -------
        numpy.ndarray
            An array containing the atom pairs that are within the specified distance threshold.
    """
    logger.info("Creating atom pairs list with a specific threshold...")
    contacts_list = []

    # Get the list of all residues.
    # A residue is identified by its name and its id.
    residues_list = list(mol.residues)

    if len(residues_list) > 1:
        # Compare residues side-by-side, between 0 and N-1 and between 1 and N
        for residue_1, residue_2 in zip(residues_list[:-1], residues_list[1:]):
            # print(f"Finding contacts for {residue_1.resid}-{residue_1.resname} and {residue_2.resid}-{residue_2.resname}")
            # Concatenate atom coordinates from both residues.
            # Example:
            #   residue_1.atoms.positions = [53.05 64.78 78.47]
            #   residue_2.atoms.positions = [54.31 65.49 78.92]
            #   ==> coords = [[53.05 64.78 78.47]
            #                 [54.31 65.49 78.92]]
            coords = np.concatenate(
                (residue_1.atoms.positions, residue_2.atoms.positions), axis=0
            )

            # Get distance between all atoms (upper-left matrix)
            # https://docs.mdanalysis.org/stable/_modules/MDAnalysis/lib/distances.html#self_distance_array
            # Result is output as a vector with (N)(N-1)/2 values.
            # Example:
            #   coords = [[53.04999924 64.77999878 78.47000122]
            #             [54.31000137 65.48999786 78.91999817]
            #             [54.         64.         74.        ]]
            #   ==> distances = [1.51466212 4.63592606 5.15000742]
            # where 1.51466212 is the distance between coords[0] and coords[1]
            #       4.63592606 is between coords[0] and coords[2]
            #       and 5.15000742 is between coords[1] and coords[2]
            distances = self_distance_array(coords)

            # Cocatenates the list of atoms ids from both residues
            atom_ids = np.concatenate((residue_1.atoms.ids, residue_2.atoms.ids))
            # Create all possible combinations between atom ids.
            pairs = np.array(list(itertools.combinations(atom_ids, 2)))
            # Create a mask for distances below a given threshold.
            # And select only atom pairs below the threshold.
            atom_pairs = pairs[distances < threshold]

            # Append atom pairs to a bigger list.
            contacts_list.append(atom_pairs)

        # Concatenate all atom pairs.
        contacts_array2 = np.concatenate(contacts_list, axis=0)

        # Remove redundant contacts
        # We have quite a lot of redundancy since we are calculating inter-contacts
        # twice for each residue: residue_1-residue_2, residue_2-residue_3, residue_3-residue_4....
        # contacts_array2 is a pair list
        contacts_array2 = np.unique(contacts_array2, axis=0)
        return contacts_array2
    else:
        # It's mean there is only one residue in this molecular system
        # So we going to return the atom pairs between the atom of this residue itself
        distances = self_distance_array(residues_list[0].atoms.positions)
        # Create all possible combinations between atom ids.
        pairs = np.array(list(itertools.combinations(residues_list[0].atoms.ids, 2)))
        # Create a mask for distances below a given threshold.
        # And select only atom pairs below the threshold.
        atom_pairs = pairs[distances < threshold]
        return atom_pairs


def get_atom_pairs_from_guess_bonds(
    molecular_system: mda.core.universe.Universe,
) -> np.ndarray:
    """This function retrieves atom pairs within a specified distance threshold from the given molecular system.
    This specified distance based on the distance between two selected atoms and their Van der Waals radius.

    References
    ----------
        https://docs.mdanalysis.org/stable/documentation_pages/core/groups.html#MDAnalysis.core.groups.AtomGroup.guess_bonds
        https://docs.mdanalysis.org/stable/documentation_pages/topology/guessers.html#MDAnalysis.topology.guessers.guess_bonds

    Parameters
    ----------
        molecular_system: mda.core.universe.Universe
            The MDAnalysis universe object representing the molecular system.

    Returns
    -------
        numpy.ndarray
            An array containing the atom pairs that are within the specified distance threshold.
    """
    logger.info(
        "Creating atom pairs list with the Van der Waals radius of each atom..."
    )

    molecular_system.atoms.guess_bonds()
    # Example:
    # bonds = <TopologyGroup containing n bonds>, where n is the number of bonds in this system
    # Each element in bonds is <Bond between: Atom n1, Atom n2>, where n1 and n2 are the indice of atom in the system
    # We can extract information of Atom n1 and Atom n2 with bond[0] and bond[1] respectively
    # bond[0] or bond[1] gives information on atom id, atom name, residue name, residue id and segid
    bonds = molecular_system.atoms.bonds
    atom_pairs = [[bond[0].id, bond[1].id] for bond in bonds]
    return np.array(atom_pairs)


def convert_atom_pairs_to_graph(
    atom_pairs: np.ndarray, mol: mda.core.universe.Universe
) -> nx.classes.graph.Graph:
    """Convert a list of pairs to a graph and its connected components.

    Reference
    ---------
    - https://stackoverflow.com/questions/4842613/merge-lists-that-share-common-elements

    Parameters
    ----------
        atom_pairs: list
            A list of pairs representing edges in the graph.
        mol: mda.core.universe.Universe
            The MDAnalysis universe object representing the molecular system.

    Returns
    -------
        networkx.classes.graph.Graph
            A graph object representing the molecular system.
    """
    logger.info("Converting atom pairs to graph...")

    graph = nx.Graph()
    # Add all atoms as single nodes.
    graph.add_nodes_from(mol.atoms.ids)

    # Add atom pairs as edges.
    graph.add_edges_from(atom_pairs)
    return graph


def add_attributes_to_nodes(
    graph: nx.classes.graph.Graph, mol_system: mda.core.universe.Universe
) -> nx.classes.graph.Graph:
    """Add molecular attributes to graph nodes.

    Attributes are taken from the molecular system (MDAnalysis universe).
    Attributes are: atom id, atom name, residue id, and residue name.

    References
    ----------
    - https://networkx.org/documentation/stable/reference/generated/networkx.classes.function.set_node_attributes.html

    Parameters
    ----------
        graph: networkx.classes.graph.Graph
            The NetworkX graph representing the molecular system.
        mol_system: MDAnalysis.core.groups.AtomGroup
            The MDAanalysis universe representing the molecular system

    Returns
    -------
        networkx.classes.graph.Graph
            The NetworkX graph representing the updated molecular system.
    """
    # logger.info(f"Adding attributes to {graph.number_of_nodes():,} nodes...")
    # logger.opt(lazy=True).debug("10 first nodes:")
    # for node_id, node_attr in sorted(graph.nodes.items())[:10]:
    #     logger.opt(lazy=True).debug(f"Node id: {node_id}")
    #     logger.opt(lazy=True).debug(f"attributes: {node_attr}")

    # Define attributes in batch: one attribute for all nodes at once.
    # This is possible because the order of nodes in the NetworkX graph
    # is the same as the order of atoms in the MDAnalysis universe.
    # Examples for note attributes after the graph is updated:
    # {'atom_id': 47, 'atom_name': 'N', 'residue_id': 3, 'residue_name': 'ALA'}
    # {'atom_id': 49, 'atom_name': 'CA', 'residue_id': 3, 'residue_name': 'ALA'}
    # {'atom_id': 51, 'atom_name': 'CB', 'residue_id': 3, 'residue_name': 'ALA'}
    nx.set_node_attributes(
        graph, dict(zip(graph.nodes, mol_system.atoms.ids)), "atom_id"
    )
    nx.set_node_attributes(
        graph, dict(zip(graph.nodes, mol_system.atoms.names)), "atom_name"
    )
    nx.set_node_attributes(
        graph, dict(zip(graph.nodes, mol_system.atoms.resids)), "residue_id"
    )
    nx.set_node_attributes(
        graph, dict(zip(graph.nodes, mol_system.atoms.resnames)), "residue_name"
    )
    # logger.opt(lazy=True).debug("10 first nodes with updated attributes:")
    # for node_id, node_attr in sorted(graph.nodes.items())[:10]:
    #     logger.opt(lazy=True).debug(f"Node id: {node_id}")
    #     logger.opt(lazy=True).debug(f"attributes: {node_attr}")
    return graph


def get_graph_components(graph: nx.classes.graph.Graph) -> list[nx.classes.graph.Graph]:
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


def get_graph_fingerprint(
    graph: nx.classes.graph.Graph,
) -> tuple[int, int, dict[str, int], list[str], dict[int, int]]:
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
                - int: The number of nodes in the graph.
                - int: The number of edges in the graph.
                - dict[str, int]: A dictionary containing the counts of each atom name present in the graph.
                - list[str]: A list containing the names of the unique residues present in the graph.
                - dict[int, int]: A dictionary containing the counts of each degree present in the graph.
    """
    nodes = graph.number_of_nodes()
    edges = graph.number_of_edges()

    atom_names = Counter(nx.get_node_attributes(graph, "atom_name").values())
    # Use sorted so the atom_names will be in the same order if it's the same molecule.
    # Otherwise it's different even if it's the same molecule
    # And the dictionnary comparaison will say it's different molecule
    atom_names = dict(sorted(atom_names.most_common()))

    # Get all residue ids and resides names pairs.
    residue_pairs = zip(
        nx.get_node_attributes(graph, "residue_id").values(),
        nx.get_node_attributes(graph, "residue_name").values(),
    )
    # Convert to dictionnary to have only one residue id (key is unique in dict).
    residue_pairs_dict = dict(residue_pairs)
    # Then extract residue names ordered by residue ids:
    residue_names = [residue_pairs_dict[key] for key in sorted(residue_pairs_dict)]

    # Exemple :
    # graph.degree = [(1, 1), (2, 3), (3, 1), (4, 2), (5, 1)]
    # dict(graph.degree) = { 1: 1, 2: 2, 3: 1, 5: 1 }
    graph_degrees_dict = dict(
        Counter([degree for _, degree in graph.degree]).most_common()
    )

    return (nodes, edges, atom_names, residue_names, graph_degrees_dict)


def get_graph_fingerprint_str(
    graph: nx.classes.graph.Graph,
    check_connectivity: bool,
) -> tuple[int, str, list[str]]:
    """Collect the tuple return by get_graph_fingerprint, to only extract
    the number of nodes, the dictionary of atom_name (that we going to convert to str so it's haschable)
    and the list of res_names.

    Parameters
    ----------
        graph: networkx.classes.graph.Graph
            The graph for which the fingerprint is to be generated.
        check_connectivity: bool
            Either we want the degree in the fingerprint or not. By default: false.

    Returns
    -------
        tuple[int, str, list[str]]:
            A tuple containing the concatenated fingerprint of the graph, which includes the following elements:
                - int: The number of nodes in the graph.
                - str: A string containing the counts of each atom name present in the graph. Exemple : "{'CA': 5, 'N': 3}"
                - list[str]: A list containing the names of the unique residues present in the graph.
    """
    (nodes, edges, atom_names, res_names, graph_degrees_dict) = get_graph_fingerprint(
        graph
    )
    if check_connectivity:
        return (nodes, edges, str(atom_names), res_names, str(graph_degrees_dict))
    else:
        return (nodes, str(atom_names), res_names)


def print_graph_fingerprint(graph: nx.classes.graph.Graph, index_graph: int):
    """Print a graph fingerprint.

    Parameters
    ----------
        graph : networkx.classes.graph.Graph
            A NetworkX graph object.
    """
    logger.debug("print groupby ... ")
    fingerprint = get_graph_fingerprint(graph)
    logger.debug(f"Graph {index_graph} fingerprint-----------------")
    logger.debug(f"- Number of nodes: {fingerprint[0]}")
    logger.debug(f"- Number of edges: {fingerprint[1]}")
    logger.debug(f"- Sorted atom names (first 50 char.): {fingerprint[2]}")
    logger.debug(f"- Sorted set of residue names: {fingerprint[3]}")
    logger.debug(f"- Node degrees dist: {fingerprint[4]}")


def get_intervals(seq: list[int]) -> list[str]:
    """Generate a list of intervals from a sorted sequence of integers.

    Resources
    ---------
        https://codereview.stackexchange.com/questions/220072/construct-intervals-from-a-sequence-of-numbers-and-vice-versa

    Parameter
    ---------
        seq: list[int])
            A sorted list of integers.

    Returns
    -------
        list[str]
            A list of strings representing intervals in the input sequence.

    Example:
        get_intervals([1, 2, 3, 6, 7, 8, 10])
        ['1-3', '6-8', '10']
    """
    starts = [x for x in seq if x - 1 not in seq]
    ends = [y for y in seq if y + 1 not in seq]
    return [str(a) + "-" + str(b) if a != b else str(a) for a, b in zip(starts, ends)]


def extract_interval(graph: nx.classes.graph.Graph) -> dict[str, list[int]]:
    """Extract residue and atom intervals from a graph and return them as a dictionary.

    Parameter
    ---------
        graph: nx.classes.graph.Graph
            A NetworkX graph with nodes that have "residue_id", "residue_name", and "atom_id" attributes.

    Returns
    -------
        dict[str, list[int]]
            A dictionary containing lists of residue IDs, residue ID intervals,
                    atom IDs, and atom ID intervals.
    """
    residue_pairs = zip(
        nx.get_node_attributes(graph, "residue_id").values(),
        nx.get_node_attributes(graph, "residue_name").values(),
    )
    residue_pairs_dict = dict(residue_pairs)
    residue_id = [key for key in sorted(residue_pairs_dict)]

    res_id = residue_id
    res_id_interval = get_intervals(residue_id)

    atom_id = nx.get_node_attributes(graph, "atom_id").values()
    atom_id_interval = get_intervals(atom_id)

    dict_res = {
        "res_id": res_id,
        "res_id_interval": res_id_interval,
        "atom_id": atom_id,
        "atom_id_interval": atom_id_interval,
    }
    return dict_res


def get_formula_based_atom_name(atom_name_dict: dict[str, int]) -> str:
    """Generates a molecular formula string, by apply guess_atom_type() from MDA on each atom names and save the occurrence.

    Ressources
    ----------
    https://docs.mdanalysis.org/2.0.0/documentation_pages/topology/guessers.html#guessing-elements-from-atom-names
    https://docs.mdanalysis.org/2.0.0/documentation_pages/topology/guessers.html#MDAnalysis.topology.guessers.guess_atom_type

    Parameters
    ----------
        atom_name_dict: dict[str, int]
            A dictionary with atom names as keys and their respective counts as values.

    Returns
    -------
        str
            A string representing the chemical formula, without hydrogens
    """
    atom_names = []
    for atom_name in atom_name_dict.keys():
        # Example:
        # mda.topology.guessers.guess_atom_type("CA") = 'C
        # mda.topology.guessers.guess_atom_type("HB1") = 'H
        atom_name = f"{mda.topology.guessers.guess_atom_type(atom_name)}"
        # print(mda.topology.guessers.guess_atom_element(atom.name))
        if atom_name != "H":
            atom_names.append(atom_name)

    atom_names = Counter(atom_names)
    sorted_atom_counts = sorted(atom_names.items())
    formula = "".join(
        f"{atom}{count}" if count > 1 else atom for atom, count in sorted_atom_counts
    )

    return formula


def count_molecule(
    graph_list: list[nx.classes.graph.Graph],
    check_connectivity: bool,
) -> dict[nx.classes.graph.Graph, dict[str, int]]:
    """Count the occurrence of molecules in a list of graphs based on their fingerprints.

    This function takes a list of graphs and counts the occurrence of each unique molecule
    based on their fingerprints, which are calculated using the get_graph_fingerprint_str function.
    
    Reference 
    ---------
    - https://stackoverflow.com/questions/46999771/comparing-a-large-number-of-graphs-for-isomorphism

    Parameters
    ----------
        graph_list: list
            A list of graph objects.
        check_connectivity: bool
            Either we want the degree in the fingerprint or not. By default: false.

    Returns
    -------
        dict
            A dictionary where keys are unique molecules (graph objects) and values are\
                their respective occurrence counts in the input list of graphs.
    """
    logger.info("Counting molecules...")
    dict_count = {}

    sorted_graphs = sorted(
        graph_list, key=lambda x: get_graph_fingerprint_str(x, check_connectivity)
    )

    for fingerprint, graph in groupby(
        sorted_graphs, key=lambda x: get_graph_fingerprint_str(x, check_connectivity)
    ):
        # fingerprint : (nb_node, nb_edge, atom_name, resname, degree)
        # graph : objet itertools that group all graph with the same fingerprint

        # A list that contain all graph with the same fingerprint
        similar_graphs = list(graph)
        nb_graph = len(similar_graphs)  # Number of graph for this fingerprint

        atom_id, res_id = [], []
        for graph in similar_graphs:
            residue_pairs = zip(
                nx.get_node_attributes(graph, "residue_id").values(),
                nx.get_node_attributes(graph, "residue_name").values(),
            )
            residue_pairs_dict = dict(residue_pairs)

            res_id.extend([key for key in sorted(residue_pairs_dict)])
            atom_id.extend(nx.get_node_attributes(graph, "atom_id").values())

        res_id_interval = get_intervals(res_id)
        atom_id_interval = get_intervals(sorted(atom_id))

        atom_names = Counter(nx.get_node_attributes(graph, "atom_name").values())
        atom_names = dict(sorted(atom_names.most_common()))
        formula = get_formula_based_atom_name(atom_names)

        # If for this fingerprint, there is only one graph
        if nb_graph == 1:
            dict_count[similar_graphs[0]] = {
                "res_id": res_id,
                "res_id_interval": res_id_interval,
                "atom_id": atom_id,
                "atom_id_interval": atom_id_interval,
                "formula_no_h": formula,
                "graph": nb_graph,
            }
        else:
            # If for this fingerprint, all the graph only have one node
            if fingerprint[0] == 1:
                dict_count[similar_graphs[0]] = {
                    "res_id": res_id,
                    "res_id_interval": res_id_interval,
                    "atom_id": atom_id,
                    "atom_id_interval": atom_id_interval,
                    "formula_no_h": formula,
                    "graph": nb_graph,
                }
            else:
                dict_count[similar_graphs[0]] = {
                    "res_id": res_id,
                    "res_id_interval": res_id_interval,
                    "atom_id": atom_id,
                    "atom_id_interval": atom_id_interval,
                    "formula_no_h": formula,
                    "graph": nb_graph,
                }
    return dict_count


def print_graph_inventory(graph_dict: dict):
    """Print graph inventory.

    Parameters
    ----------
        graph_dict: dict
            A dictionary with graphs as keys and counts (numbers of graphs) as values.
    """
    logger.info("Molecular inventory:")
    total_molecules_count = 0
    for graph_idx, (graph, key) in enumerate(graph_dict.items(), start=1):
        logger.info(f"Molecule {graph_idx:,} ----------------")

        if len(key) == 6:
            (_, res_id_interval, _, _, _, count) = key.values()
        else:
            (_, res_id_interval, _, _, name, count, _) = key.values()
            logger.info(f"- name: {name}")
        # else:
        #     (_, res_id_interval, _, _, name, count, _, _, _) = key.values()
        #     logger.info(f"- name: {name}")

        logger.info(f"- number of atoms: {graph.number_of_nodes():,}")
        logger.info(f"- number of molecules: {count:,}")

        logger.debug("- residue ids:")
        for i in range(min(10, len(res_id_interval))):
            logger.debug(f"\t({res_id_interval[i][:20]})")

        res_names = set(sorted(nx.get_node_attributes(graph, "residue_name").values()))
        logger.debug(f"- res names: {res_names}")

        total_molecules_count += count
    logger.success(f"{total_molecules_count:,} molecules in total")


def print_first_atoms(
    mda_universe: mda.core.universe.Universe, number_of_atoms: int = 10
):
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
        logger.debug(
            f"Res. id: {atom.residue.resid} | res. name: {atom.residue.resname} | "
            f"at. name: {atom.name} | at. id: {atom.id} | "
            f"at. coord.: {string_coords}"
        )


def read_structure_file_remove_hydrogens(file_path: str) -> mda.core.universe.Universe:
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
    molecule_without_h.write(without_h_file_path, reindex=False)
    return molecule_without_h


def remove_hydrogene(filename: str) -> mda.core.universe.Universe:
    """Removes hydrogen atoms from a molecular system.

    Parameters
    ----------
        filename : str
            Path to the file containing the molecular structure.

    Returns
    -------
        mda.core.universe.Universe
            MDAnalysis Universe object representing the molecular system without hydrogen atoms.
    """
    molecule = mda.Universe(filename)
    logger.info(f"Found {len(molecule.atoms):,} atoms in {filename}")

    # Remove hydrogene from the system
    logger.info("Removing hydrogen atoms...")
    mol = molecule.select_atoms("not (name H* or name [123456789]H*)")
    filename_tmp = f"{Path(filename).stem}_without_H{Path(filename).suffix}"
    # Write the new system in a new file
    mol.write(filename_tmp, reindex=False)
    logger.debug(f" New structure file without hydrogens : {filename_tmp}")

    # We need to read structure from disk to be extra sure hydrogen atoms are removed.
    mol = mda.Universe(filename_tmp)
    logger.info(f"{len(mol.atoms):,} atoms remaining")
    return mol


def find_ion_solvant(
    molecule: dict,
    universe: mda.core.universe.Universe,
    counts: dict,
    solvant_or_ion: str,
) -> tuple[mda.core.universe.Universe, dict[nx.classes.graph.Graph, dict[str, int]]]:
    """Counts and removes ions or solvents from the MDAnalysis Universe.

    Parameters
    ----------
        molecule : dict
            Dictionary containing information about the atoms (ion or solvant) to be removed (name, res_name and atom_name).
        universe : MDAnalysis.core.universe.Universe
            MDAnalysis Universe object representing the system.
        counts : dict
            Dictionary to store the counts of ions or solvents.
        solvant_or_ion : str
            Information about the type of the molecule

    Returns
    -------
        MDAnalysis.core.universe.Universe
            MDAnalysis Universe object containing only non-ion or non-solvent atoms.
        dict
            Dictionary containing the counts of removed atoms.
            With the graph (representing this molecule) as key,
            and in the value :
                - the atom_id of the first atom (for each molecule)
                - the atom_id of the last atom (for each molecule)
                - the name of this molecule (collected from the dictionary in mol_def.py)
                - the occurence if this molecule in this system
                - boolean key for ion, solvant and lipid
    """
    (name, res_name, atom_names) = molecule.values()

    # To select the ion (or solvant) by their res_name and all their atom_name (if there are multiple)
    selection = f"resname {res_name} and (name {' or name '.join(atom_names)})"
    selected_atoms = universe.select_atoms(selection)

    # Collect all resids from each residues selected, to remove it from the univers
    selected_res_ids = [str(residue.resid) for residue in selected_atoms.residues]
    res_count = len(selected_res_ids)

    # Only change the Universe and the count dictionnary if this res_name is in this Universe
    # Otherwise we retrun the same Universe and same count dictionnary
    if res_count > 0:
        list_graph = []
        # For each residues in the selection, create a graph in the same format as a molecule (but added to that a key 'name')
        # which give us direct acces of their composition
        for index_resID in selected_atoms.residues:
            graph = nx.Graph()
            graph.add_nodes_from(index_resID.atoms.ids)
            graph = add_attributes_to_nodes(graph, index_resID)
            list_graph.append(graph)

        atom_id, res_id = [], []
        for subgraph in list_graph:
            residue_pairs = zip(
                nx.get_node_attributes(subgraph, "residue_id").values(),
                nx.get_node_attributes(subgraph, "residue_name").values(),
            )
            residue_pairs_dict = dict(residue_pairs)

            res_id.extend([key for key in sorted(residue_pairs_dict)])
            atom_id.extend(nx.get_node_attributes(subgraph, "atom_id").values())

        res_id_interval = get_intervals(res_id)
        atom_id_interval = get_intervals(atom_id)

        if solvant_or_ion == "ion":
            molecular_type = "ion"
        else:
            molecular_type = "solvant"

        counts[list_graph[0]] = {
            "res_id": res_id,
            "res_id_interval": res_id_interval,
            "atom_id": atom_id,
            "atom_id_interval": atom_id_interval,
            "name": name,
            "graph": res_count,
            "molecular_type": molecular_type,
        }

        # Here we remove all the resIDS (from selected_res_ids) from this universe
        for interval in res_id_interval:
            start_end = interval.split("-")
            if len(start_end) == 1:
                selection = f"not (resname {res_name} and (name {' or name '.join(atom_names)}) and resid {start_end[0]})"
            else:
                selection = f"not (resname {res_name} and (name {' or name '.join(atom_names)}) and resid {start_end[0]}:{start_end[1]})"

            universe = universe.select_atoms(f"{selection}")
    return (universe, counts)


def count_remove_ion_solvant(
    universe: mda.core.universe.Universe,
    input_filepath: str,
) -> tuple[mda.core.universe.Universe, dict[nx.classes.graph.Graph, dict[str, int | str]]]:
    """Count and remove ions and solvents from the MDAnalysis Universe return by
    the function find_ion_solvant().

    Parameters
    ----------
        universe : MDAnalysis.core.universe.Universe
            MDAnalysis Universe object representing the system.
        input_filepath : str
            Path to the input file.

    Returns
    -------
        tuple
            Containing :
                - the new Universe without ions and solvants.
                - a dictionary where
                    - the key is a graph
                    - and the value is an other dictionary with: atom_start, atom_end, name of the ion-solvant, the counts of removed ions-solvant
    """
    counts = {}

    logger.info("Searching ions...")
    for ion in mol_def.IONS_LIST:
        universe, counts = find_ion_solvant(ion, universe, counts, "ion")

    logger.info("Searching solvant molecules...")
    for solvant in mol_def.SOLVANTS_LIST:
        # Solvant methanol and the residue methionine could be confused
        # because their share the same residue name 'MET'.
        # We ignore methanol for now.
        if solvant["res_name"] == "MET":
            continue
        universe, counts = find_ion_solvant(solvant, universe, counts, "solvant")

    # Write the new universe without ions and solvant into a new file
    output_file = f"{Path(input_filepath).stem}_without_H_ions_solvant{Path(input_filepath).suffix}"
    universe.atoms.write(output_file, reindex=False)
    logger.debug(
        f" New structure file without hydrogens, ions and solvants : {output_file}"
    )

    universe_clean = mda.Universe(output_file)

    # Print which ion and solvant we find, and how many
    for molecule, dict_count in counts.items():
        name = dict_count.get("name")
        count = dict_count.get("graph")
        res_name = " ".join(
            set(nx.get_node_attributes(molecule, "residue_name").values())
        )
        logger.success(f"Found: {count} {name} ({res_name})")

    # Check if there is other residue with the resname SOL in the updated MDAnalysis.core.universe.Universe
    selected_atoms = universe_clean.select_atoms("resname SOL")
    count = len(selected_atoms.residues)
    logger.info(f"{count} residues SOL remaining")

    logger.info(f"{len(universe_clean.atoms):,} atoms remaining")
    return (universe_clean, counts)


def check_overlapping_residue_between_graphs(graph_list: list[nx.classes.graph.Graph]):
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
        # Here it only compare the residue id
        # But in some case, the residue id is reinitialize
        # So some molecule will have the same id but it's not the same
        # res_id_set = set((nx.get_node_attributes(graph, "residue_id").values()))

        # Here I add the residue name to the comparaison
        # So we see the overlapping with the residue id and the residue name
        res_id_set = set((nx.get_node_attributes(graph, "residue_id").values()))
        res_id_name_list = [(nx.get_node_attributes(graph, "residue_name").values())]
        res_id_set = set((tuple(res_id_set), tuple(res_id_name_list)))
        # print(res_id_set)

        res_id_intersect = res_id_set_all.intersection(res_id_set)
        res_id_set_all.update(res_id_set)
        if res_id_intersect:
            for res_id in res_id_intersect:
                logger.critical(f"Residue id {res_id} is found in multiple graphs")
                res_id_common.append(res_id)

    if not res_id_common:
        logger.success("No overlapping residue found")
        return None
    else:
        res_id_common_set = set(res_id_common)
        # for res_id in res_id_common_set:
        #     for graph_id, graph in enumerate(graph_list):
        #         res_id_set = set((nx.get_node_attributes(graph, "residue_id").values()))
        #         if res_id in res_id_set:
        #             logger.error(f"Residue id {res_id} is found in graph id {graph_id}")
        return res_id_common_set
        # raise Exception("Some residue id are found in multiple graphs")


def is_protein(graph: nx.classes.graph.Graph) -> bool:
    """Check if the molecule represented by the graph is a protein.

    This function checks whether the graph represents a protein molecule by
    verifying if at least 3 residues names are in the the list of residues names given by MDAnalysis.

    Parameters
    ----------
        graph: networkx.Graph
            The input graph representing the molecule.

    Returns
    -------
        bool
            True if the molecule is a protein, False otherwise.
    """
    set_key_amino_acid_mda = set(mol_def.AMINO_ACID_DICT.keys())
    set_res_name_graph = set(nx.get_node_attributes(graph, "residue_name").values())
    return len(set_key_amino_acid_mda.intersection(set_res_name_graph)) > 3


def extract_protein_sequence(graph: nx.classes.graph.Graph) -> dict[str, int]:
    """Extract the protein sequence from a graph.

    This function extracts the protein sequence from the molecule represented
    by the input graph. By getting all the residue name sorted by their residue ids (to be sure it's in the right order).

    Parameters
    ----------
        graph: networkx.Graph
            The input graph representing the molecule.

    Returns
    -------
        dict
            A dictionary containing the following keys:
                - 'sequence': str
                    The protein sequence extracted from the molecule.
                - 'nb_res': int
                    The number of residue in this protein.
    """
    logger.info("Extracting protein sequence...")
    info_seq = {}
    protein_sequence = []

    residue_pairs = zip(
        nx.get_node_attributes(graph, "residue_id").values(),
        nx.get_node_attributes(graph, "residue_name").values(),
    )
    residue_pairs_dict = dict(residue_pairs)
    residue_names = [residue_pairs_dict[key] for key in sorted(residue_pairs_dict)]
    for resname in residue_names:
        protein_sequence.append(mol_def.AMINO_ACID_DICT.get(resname, "?"))

    info_seq["sequence"] = "".join(protein_sequence)
    info_seq["nb_res"] = len(protein_sequence)
    return info_seq


def export_protein_sequence_into_FASTA(
    protein_sequence_dict: dict[int, dict[str, int]], filepath_name: str
):
    """Export the protein sequences into a FASTA file.

    Parameters
    ----------
        protein_sequence_dict : dict
            A dictionary containing the protein sequences, where the keys are
            identifiers and the values are dictionaries with the following keys:
                - 'sequence': str
                    The protein sequence.
                - 'nb_res': int
                    The number of residues in the protein sequence.
        filepath_name : str
            The filepath for the output FASTA file.
    """
    logger.info("Converting into FASTA file...")
    with open(filepath_name, "w") as file:
        for info_seq in protein_sequence_dict.values():
            seq, nb_res = info_seq.values()
            # For only have 80 residues for each line
            seq = [seq[i : i + 80] for i in range(0, len(seq), 80)]
            content = f">Protein: {nb_res} residues\n" + "\n".join(seq)
            file.write(f"{content}\n")
    logger.debug(f"FASTA filename : {filepath_name}")


def is_lipid(
    resolution: str, graph: nx.classes.graph.Graph, dict_count: dict[str, int|str]
) -> bool:
    """Determines if the given graph represents a lipid.

    Parameters
    ----------
        resolution: str
            The resolution of the molecular system.
        graph: nx.classes.graph.Graph
            The molecular graph with nodes containing attributes, including "residue_name".
        dict_count: dict[str, str]
            A dictionary containing molecular information with keys such as "formula_no_h".

    Returns
    -------
        bool
            True if the graph represents a lipid according to the predefined database, False otherwise.
    """
    res_name_graph = set(nx.get_node_attributes(graph, "residue_name").values())
    res_name_graph = res_name_graph.pop()

    if resolution == "AA":
        lipid_csml_charmm_gui = CSML_CHARMM_GUI[
            CSML_CHARMM_GUI["Category"].str.contains("lipid", case=False, na=False)
        ]
        if "formula_no_h" in dict_count.keys():
            formula_graph = dict_count["formula_no_h"]
            selected_row = lipid_csml_charmm_gui.loc[
                (lipid_csml_charmm_gui["Alias"].apply(lambda x: res_name_graph in x)) 
                & (lipid_csml_charmm_gui["Formula"] == formula_graph)]

            # Check if the selection match a row with this condition
            if not selected_row.empty:
                dict_count["name"] = list(selected_row["Name"].values)
                return True
    else:
        lipid_MAD = MAD_DB[
            MAD_DB["Category"].str.contains("Lipids", case=False, na=False)
        ]
        selected_row = lipid_MAD.loc[(lipid_MAD["Alias"] == res_name_graph)]
        if not selected_row.empty:
            dict_count["name"] = list(selected_row["Name"].values)
            return True
    return False


def guess_resolution(
    molecular_system: mda.core.universe.Universe,
    number_of_res: int = 5,
    threshold: float = 2.0,
) -> str:
    """Finds the resolution of the molecular system.

    - Select 3 residues with more than 2 atoms.
    - Compute distance between atoms of each residue.
    - System is all-atom ("AA") if the minimum distance between atoms of each residue is less than 2 Å.
    - System is coarse-grained ("CG") otherwise.

    Parameters
    ----------
        molecular_system: MDAnalysis.Universe
            The molecular universe without hydrogen ions and solvent.
        number_of_res: int
            The number of residues to collect. By default: 5.
        threshold: float
            Threshold distance to consider the system as all-atom. By default: 2.0 Å.

    Returns
    -------
        str
            The resolution of the molecular system:
            either "AA" (all_atom) or "CG" (coarse-grain).
    """
    residues_to_analyze = [
        residue
        for residue in molecular_system.atoms.residues
        if len(residue.atoms) >= 2
    ]
    # System is all-atom ("AA") is at least one residue
    # has inter-atomic distance below threshold.
    for residue in residues_to_analyze[:number_of_res]:
        distances = self_distance_array(residue.atoms.positions)
        if sum(distances < threshold):
            return "AA"
    return "CG"


def get_git_last_commit_date() -> str:
    """Get the last commit date from the git repository."""
    try:
        command = "git show --no-patch --no-notes --pretty='%cI' HEAD"
        git_date = subprocess.check_output(command.split()).decode("ascii").strip()
        return git_date
    except subprocess.CalledProcessError:
        return ""


def get_git_last_commit_hash() -> str:
    """Get the last commit hash from the git repository."""
    try:
        command = "git show --no-patch --no-notes --pretty='%H' HEAD"
        git_hash = subprocess.check_output(command.split()).decode("ascii").strip()
        return git_hash
    except subprocess.CalledProcessError:
        return ""


def export_inventory(
    graph_count_dict: dict[nx.classes.graph.Graph, dict[str, int]],
    resolution: str,
    filename: str,
    execution_time: float,
    overlap_residue: None | set[int],
) -> str:
    """Exports inventory data from a dictionary of graph objects and their associated information about molecule into JSON file.

    Parameters:
    ----------
        graph_count_dict : dict
            A dictionary where each key is a NetworkX graph object, and each value is another dictionary containing
            various attributes and metadata about the graph, including whether it represents a protein, lipid, ion,
            or solvent.
        resolution: str
            The resolution of the molecular system.
        filename: str
            Path to the file containing the molecular structure.
    """
    list_dict_molecule = []
    for index_graph, (graph, information) in enumerate(
        graph_count_dict.items(), start=1
    ):
        formula, molecular_type, protein_sequence, remark_message, comment = (
            "",
            "unknown",
            "",
            "",
            "",
        )
        putative_pdb, putative_name = [], []

        residue_pairs = zip(
            nx.get_node_attributes(graph, "residue_id").values(),
            nx.get_node_attributes(graph, "residue_name").values(),
        )
        residue_pairs_dict = dict(residue_pairs)
        residue_names = [residue_pairs_dict[key] for key in sorted(residue_pairs_dict)]
        residue_names = " ".join(residue_names)

        atom_names = Counter(nx.get_node_attributes(graph, "atom_name").values())
        atom_names = dict(sorted(atom_names.most_common()))

        if "formula_no_h" in information.keys():
            formula = information["formula_no_h"]
        if "molecular_type" in information.keys():
            molecular_type = information["molecular_type"]
            if molecular_type == "protein":
                protein_sequence = information["protein_sequence"]
                if "putative_pdb" in information.keys():
                    putative_pdb = information["putative_pdb"]
            elif molecular_type in ["lipid", "ion", "solvant"]:
                putative_name = information["name"]
        if "comment" in information.keys():
            comment = information["comment"]

        dict_inventory = {
            "id": index_graph,
            "number_of_atoms": graph.number_of_nodes(),
            "number_of_molecules": information["graph"],
            "residue_names": residue_names,
            "residue_ids": " ".join(information["res_id_interval"]),
            "formula_without_h": formula,
            "molecular_type": molecular_type,
            "protein_sequence": protein_sequence,
            "putative_pdb": putative_pdb,
            "putative_name": putative_name,
            "comment": comment,
        }
        list_dict_molecule.append(dict_inventory)

        date_time = f"{datetime.datetime.now():%Y-%m-%d_%H-%M-%S}"

        # current_dir = Path.cwd()
        # absolute_file_path = current_dir / filename
        # relative_path = absolute_file_path.relative_to(current_dir)
        relative_path = filename

        if overlap_residue:
            remark_message = [
                f"Residue {x} has been splitted into multiple molecules. This should be wrong."
                for x in overlap_residue
            ]

        final_dict = {
            "inventory": sorted(
                list_dict_molecule, key=lambda x: x["number_of_molecules"]
            ),
            "resolution": resolution,
            "date": date_time,
            "execution_time_in_sec": f"{execution_time:.2f}",
            "remark": remark_message,
            "file_path": str(relative_path),
            "file_md5sum": hashlib.md5(open(filename, "rb").read()).hexdigest(),
            "git_last_commit_date": get_git_last_commit_date(),
            "git_last_commit_hash": get_git_last_commit_hash(),
        }

    logger.info("Exporting inventory into JSON file...")
    filename_JSON = f"{date_time}_grodecoder_{Path(filename).stem}.json"
    out_file = open(filename_JSON, "w")
    json.dump(final_dict, out_file, indent=4)
    out_file.close()
    logger.info(f"JSON filename : {filename_JSON}")
    return filename_JSON


def is_met(graph: nx.classes.graph.Graph) -> bool:
    """Determines if a given graph represents a methanol (MET).

    Parameters
    ----------
        graph: nx.classes.graph.Graph
            A NetworkX graph representing a molecular structure.

    Returns
    -------
        bool
            True if the graph represents a MET solvant, False otherwise.
    """
    res_name = (set(nx.get_node_attributes(graph, "residue_name").values())).pop()
    nb_atom = graph.number_of_nodes()

    if res_name == "MET" and nb_atom == 2:
        return True
    return False


def main(
    input_file_path: str,
    check_connectivity: bool = False,
    bond_threshold: str | float = "auto",
    query_pdb=False,
):
    """Excute the main function for analyzing a .gro file.

    Parameters
    ----------
        input_file_path: str
            Filepath of the .gro file we want to analyzed
        check_connectivity: boolean
            If we want to add degrees and the number of edges of each graph in their fingerprint. By default at False.
        bond_threshold : str | float
            Choose the method used to get the atom pairs, there is two options: 'auto' if we don't know the resolution of the system and we detect it or enter a threshold (a positiv float number) if we know it's a corse grain model.
        query_pdb: boolean
            If we want to have informations (PDB ID, name, organism) about the protein identified in the PDB API. By default at False.
    """
    start_time = time.perf_counter() 

    molecular_system = remove_hydrogene(input_file_path)
    molecular_system, count_ion_solvant = count_remove_ion_solvant(
        molecular_system,
        input_file_path,
    )

    resolution = guess_resolution(molecular_system)
    logger.info(f"Molecular resolution: {resolution}")
    if bond_threshold == "auto":
        if resolution == "AA":
            atom_pairs = get_atom_pairs_from_guess_bonds(molecular_system)
        else:
            atom_pairs = get_atom_pairs_from_threshold(molecular_system, 5.0)
    else:
        atom_pairs = get_atom_pairs_from_threshold(molecular_system, bond_threshold)

    graph_return = convert_atom_pairs_to_graph(atom_pairs, molecular_system)

    graph_with_node_attributes = add_attributes_to_nodes(graph_return, molecular_system)
    graph_list = get_graph_components(graph_with_node_attributes)

    overlap_residue = check_overlapping_residue_between_graphs(graph_list)

    graph_count_dict = count_molecule(graph_list, check_connectivity)

    for index_graph, graph in enumerate(graph_count_dict.keys(), start=1):
        print_graph_fingerprint(graph, index_graph)

    graph_count_dict.update(count_ion_solvant)
    print_graph_inventory(graph_count_dict)

    filename = Path(input_file_path).stem

    protein_sequence_dict = {}
    for index_graph, (graph, key) in enumerate(graph_count_dict.items(), start=1):
        if is_met(graph):
            graph_count_dict[graph]["molecular_type"] = "solvant"
            graph_count_dict[graph]["name"] = "organic solvant methanol/OPLS"

        elif is_protein(graph):
            sequence_nbres = extract_protein_sequence(graph)
            sequence = sequence_nbres["sequence"]

            protein_sequence_dict[index_graph] = sequence_nbres
            graph_count_dict[graph]["molecular_type"] = "protein"
            graph_count_dict[graph]["protein_sequence"] = sequence

            list_dict_info_pdb = []
            if query_pdb:
                results = search_into_PDB.API_PDB_search_based_sequence(sequence)
                if not results:
                    graph_count_dict[graph][
                        "comment"
                    ] = "No corresponding structure found in the PDB"
                for pdb_id in results:
                    list_dict_info_pdb.append(
                        search_into_PDB.get_info_one_pdb_id(pdb_id)
                    )
                graph_count_dict[graph]["putative_pdb"] = list_dict_info_pdb
        elif is_lipid(resolution, graph, key):
            graph_count_dict[graph]["molecular_type"] = "lipid"
    export_protein_sequence_into_FASTA(protein_sequence_dict, f"{filename}.fasta")

    # execution_time in seconds
    execution_time = time.perf_counter() - start_time
    logger.info(f"execution_time: {execution_time}")
    JSON_filepath = export_inventory(
        graph_count_dict, resolution, input_file_path, execution_time, overlap_residue
    )
    return JSON_filepath


def is_a_structure_file(filepath: str) -> str:
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

    if filename.suffix not in (".gro", ".pdb"):
        raise argparse.ArgumentTypeError(f"{filepath} is not a .gro or .pdb file.")
    return filepath


def is_a_valid_threshold(threshold: str) -> str | float:
    """Check if the given threshold is permitted).

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
    if threshold == "auto":
        return threshold
    try:
        threshold_as_float = float(threshold)
    except argparse.ArgumentTypeError:
        raise argparse.ArgumentTypeError("Argument should 'auto' or a number")
    if not threshold_as_float > 0.0:
        raise argparse.ArgumentTypeError("Argument should be > 0")
    return threshold_as_float


def parse_arg() -> argparse.Namespace:
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
    parser = argparse.ArgumentParser(
        prog="grodecoder",
        description="Extract molecules from a structure file (.gro, .pdb).",
        usage="grodecoder.py [-h] --input structure_file [--drawgraph]",
    )
    parser.add_argument(
        "--input",
        type=is_a_structure_file,
        help="structure file path (.gro, .pdb)",
        required=True,
    )
    parser.add_argument(
        "--checkconnectivity",
        help="Add edges and degre in the fingerprint. Default: False.",
        default=False,
        action="store_true",
    )
    parser.add_argument(
        "--bondthreshold",
        type=is_a_valid_threshold,
        help="Choose the method to calculate the atom pairs. If we know the resolution of the system is coarse grain enter a threshold (a positiv float number) or we don't know so choose 'auto'",
        default="auto",
    )
    parser.add_argument(
        "--querypdb",
        help="Add PDB id and their putative name in the JSON file for the protein. Default: False.",
        default=False,
        action="store_true",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_arg()
    main(
        args.input,
        check_connectivity=args.checkconnectivity,
        bond_threshold=args.bondthreshold,
        query_pdb=args.querypdb,
    )
