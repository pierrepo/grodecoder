"""Extract each molecule of a GRO file and print their occurence.

Usage:
    import grodecoder as gd
"""

__authors__ = "Karine DUONG"
__contact__ = "karine.duong@etu.u-paris.fr"
__copyright__ = "IBPC"
__date__ = "2024-03-18"


from collections import Counter
import itertools
from itertools import groupby
from pathlib import Path
import re

import argparse
from loguru import logger
import MDAnalysis as mda
from MDAnalysis.analysis.distances import contact_matrix, self_distance_array
import matplotlib.pyplot as plt
import numpy as np
import networkx as nx
from networkx.algorithms.components.connected import connected_components
from scipy.sparse import triu
from scipy.spatial.distance import cdist

import mol_def
import fasta_seq_into_PDBID as fsiPDB


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
        for index_j, pos_j in enumerate(position[index_i + 1 :], start=index_j + 1):
            tmp[index_i][index_j] = cdist([pos_i], [pos_j])
    return tmp


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


def get_atom_pairs2(mol, threshold):
    contacts_list = []

    # Get the list of all residues.
    # A residue is identified by its name and its id.
    residues_list = list(mol.residues)

    # Compare residues side-by-side, between 0 and N-1 and between 1 and N
    for residue_1, residue_2 in zip(residues_list[:-1], residues_list[1:]):
        # print(f"Finding contacts for {residue_1.resid}-{residue_1.resname} and {residue_2.resid}-{residue_2.resname}")
        # Concatenate atom coordinates from both residues.
        coords = np.concatenate(
            (residue_1.atoms.positions, residue_2.atoms.positions), axis=0
        )

        # Get distance between all atoms (upper-left matrix)
        # https://docs.mdanalysis.org/stable/_modules/MDAnalysis/lib/distances.html#self_distance_array
        # Result is output as a vector with (N)(N-1)/2 values.
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


def convert_atom_pairs_to_graph(atom_pairs, mol):
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
    # graph.add_nodes_from(list(range(0, total_number_of_atoms)))
    graph.add_nodes_from(mol.atoms.ids)

    # Add atom pairs as edges.
    graph.add_edges_from(atom_pairs)
    return graph


def add_attributes_to_nodes(graph, mol_system):
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
    logger.info(f"Adding attributes to {graph.number_of_nodes():,} nodes...")
    logger.opt(lazy=True).debug("10 first nodes:")
    for node_id, node_attr in sorted(graph.nodes.items())[:10]:
        logger.opt(lazy=True).debug(f"Node id: {node_id}")
        logger.opt(lazy=True).debug(f"attributes: {node_attr}")
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
    logger.opt(lazy=True).debug("10 first nodes with updated attributes:")
    for node_id, node_attr in sorted(graph.nodes.items())[:10]:
        logger.opt(lazy=True).debug(f"Node id: {node_id}")
        logger.opt(lazy=True).debug(f"attributes: {node_attr}")
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
    # res_names = " ".join(
    #     sorted(set((nx.get_node_attributes(graph, "residue_name").values())))
    # )
    res_names = set()
    for res_name in set((nx.get_node_attributes(graph, "residue_name").values())):
        res_names.add(mol_def.AMINO_ACID_DICT.get(res_name, res_name))

    graph_degrees = Counter(dict(graph.degree).values())
    degree_dist = " ".join(
        [f"{key}:{value}" for key, value in sorted(graph_degrees.items())]
    )

    return (nodes, edges, atom_names, res_names, degree_dist)
    # return (nodes, edges, atom_names, res_names)
    # return (nodes, atom_names, res_names)


def print_graph_fingerprint(graph, index_graph):
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

        # A list that contain all graph with the same fingerprint
        similar_graphs = list(graph)
        nb_graph = len(similar_graphs)  # Number of graph for this fingerprint

        atom_start, atom_end = [], []
        for graph in similar_graphs:
            atom_start.append(min(nx.get_node_attributes(graph, "atom_id").values()))
            atom_end.append(max(nx.get_node_attributes(graph, "atom_id").values()))

        if nb_graph == 1:  # For this fingerprint, there is only one graph
            dict_count[similar_graphs[0]] = {
                "atom_start": atom_start,
                "atom_end": atom_end,
                "graph": nb_graph,
            }
        else:
            # For this fingerprint, all the graph only have one node
            if fingerprint[0] == 1:
                dict_count[similar_graphs[0]] = {
                    "atom_start": atom_start,
                    "atom_end": atom_end,
                    "graph": nb_graph,
                }
            else:
                dict_count[similar_graphs[0]] = {
                    "atom_start": atom_start,
                    "atom_end": atom_end,
                    "graph": nb_graph,
                }
    return dict_count


def get_ion_solvant_info(residue_name, target_field):
    """Get information about the specified molecule from the list of ions or solvants.

    Parameters
    ----------
        residue_name: str
            The name of the molecule to search for.
        target_field: str
            The information to retrieve for the molecule, such as 'name' or 'atom_names'.

    Returns
    -------
        str
            The value of the specified `reseach_term` for the molecule, if found.
    """
    for molecule in mol_def.IONS_LIST + mol_def.SOLVANTS_LIST:
        if molecule["res_name"] == residue_name:
            return molecule[target_field]
    return None


def print_graph_inventory(graph_dict):
    """Print graph inventory.

    Parameters
    ----------
        graph_dict: dict
            A dictionary with graphs as keys and counts (numbers of graphs) as values.
    """
    logger.info("Molecules inventory:")
    total_molecules_count = 0
    for graph_idx, (graph, key) in enumerate(graph_dict.items(), start=1):
        logger.info(f"Molecule {graph_idx:,} ----------------")

        if isinstance(graph, str):
            count = key
            logger.info(
                f"- number of atoms: {len(get_ion_solvant_info(graph, 'atom_names'))}"
            )
            logger.info(f"- number of molecules: {count:,}")
            logger.debug(f"- 20 first atom names: {graph}")
            logger.debug(f"- res names: {graph}")
        elif isinstance(graph, nx.Graph):
            atom_start, atom_end, count = key.values()
            logger.info(f"- number of atoms: {graph.number_of_nodes():,}")
            logger.info(f"- number of molecules: {count:,}")

            # logger.info(f"- tuple of min and max atom_id for each graph:")
            # for i in range(min(20, len(atom_start))):
            #     logger.info(f"\t({atom_start[i]:,} -- {atom_end[i]:,})")

            atom_names = list(
                sorted(nx.get_node_attributes(graph, "atom_name").values())
            )
            atom_names_str = " ".join(atom_names[:20])
            logger.debug(f"- 20 first atom names: {atom_names_str}")

            res_names = set(
                sorted(nx.get_node_attributes(graph, "residue_name").values())
            )
            logger.debug(f"- res names: {res_names}")
        else:
            raise Exception("graph is neither str or nx.Graph")

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
            # To replace number by empty space, for only keep atom name
            atom_name = re.sub(r"\d", "", graph.nodes[node]["atom_name"])
            if atom_name in ("C", "CA", "CB", "CD", "CE", "CG", "CH", "CZ"):
                node_colors.append("black")
            elif atom_name in ("O", "OD", "OE", "OG", "OH", "OT", "OW"):
                node_colors.append("red")
            elif atom_name in ("N", "ND", "NE", "NH", "NZ"):
                node_colors.append("blue")
            else:
                node_colors.append("green")  # Default color for other labels
        nx.draw(
            graph,
            node_color=node_colors,
            node_size=75,
            with_labels=True,
            labels=nx.get_node_attributes(graph, "atom_name"),
            edge_color="grey",
        )
    else:
        nx.draw(
            graph,
            node_color="green",
            node_size=75,
            with_labels=True,
            labels=nx.get_node_attributes(graph, "atom_name"),
            edge_color="grey",
        )
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
        logger.debug(
            f"Res. id: {atom.residue.resid} | res. name: {atom.residue.resname} | "
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
    molecule_without_h.write(without_h_file_path, reindex=False)
    return molecule_without_h


def remove_hydrogene(filename):
    molecule = mda.Universe(filename)
    logger.info(f"Found {len(molecule.atoms):,} atoms")

    # Remove hydrogene from the system
    logger.info("Removing H atoms...")
    mol = molecule.select_atoms("not (name H*)")
    filename_tmp = f"{Path(filename).stem}_without_H{Path(filename).suffix}"
    # Write the new system in a new file
    mol.write(filename_tmp, reindex=False)

    # We need to read structure from disk to be extra sure hydrogen atoms are removed.
    mol = mda.Universe(filename_tmp)
    logger.info(f"Found {len(mol.atoms):,} atoms remaining")
    return mol


def find_ion_solvant(atoms, universe, counts):
    """Counts and removes ions or solvents from the MDAnalysis Universe.

    Parameters
    ----------
        atoms : dict
            Dictionary containing information about the atoms (ion or solvant) to be removed (res_name and atom_name).
        universe : MDAnalysis.core.universe.Universe
            MDAnalysis Universe object representing the system.
        counts : dict
            Dictionary to store the counts of ions or solvents.

    Returns
    -------
        dict
            Dictionary containing the counts of removed atoms.
        MDAnalysis.core.universe.Universe
            MDAnalysis Universe object containing only non-ion or non-solvent atoms.
    """
    res_name = atoms["res_name"]
    atom_names = atoms["atom_names"]

    # To select the ion (or solvant) by their res_name and all their atom_name (if there are multiple)
    selection = f"resname {res_name} and (name {' or name '.join(atom_names)})"

    selected_atoms = universe.select_atoms(selection)
    # Because methanol and methionine can get confuse with their res_name "MET"
    # So if the residue in this selection don't have 6 atoms (or 2 if we remove the hydrogen), it's not a methanol
    # So we remove this residue from the selection
    list_resid_methionine = set()
    if res_name == "MET":
        for index_res in selected_atoms.residues:
            if len(index_res.atoms) != 2:
                list_resid_methionine.add(index_res.resid)
        if len(list_resid_methionine) != 0:
            tmp_select = f"{selection} and not resid {' and not resid '.join(map(str, list_resid_methionine))}"
            selected_atoms = universe.select_atoms(tmp_select)

    count = len(selected_atoms.residues)

    if count > 0:
        counts[res_name] = count
        # Remove these ions or solvents from the universe
        universe = universe.select_atoms(f"not ({selection})")
    return (universe, counts)


def count_remove_ion_solvant(universe, input_filepath):
    """Count and remove ions and solvents from the MDAnalysis Universe return by
    the function count_remove_ion_solvant_aux().

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
                - the counts of removed ions
                - solvents, and the new Universe without ions and solvents.
    """
    counts = {}

    logger.info("Searching ions...")
    for ion in mol_def.IONS_LIST:
        universe, counts = find_ion_solvant(ion, universe, counts)

    logger.info("Searching solvant molecules...")
    for solvant in mol_def.SOLVANTS_LIST:
        universe, counts = find_ion_solvant(solvant, universe, counts)

    # Write the new universe without ions and solvant into a new file
    output_file = f"{Path(input_filepath).stem}_without_H_ions_solvant{Path(input_filepath).suffix}"
    universe.atoms.write(output_file, reindex=False)
    for molecule, count in counts.items():
        res_name = get_ion_solvant_info(molecule, "name")
        logger.success(f"Found: {count} {res_name} ({molecule})")

    universe_clean = mda.Universe(output_file)
    selected_atoms = universe.select_atoms("resname SOL")
    count = len(selected_atoms.residues)
    logger.info(f"{count} residues SOL remaining")

    logger.info(f"Found {len(universe_clean.atoms):,} atoms remaining")
    return (universe_clean, counts)


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
        res_id_set_all.update(res_id_set)
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
        atom_names = " ".join(
            sorted(nx.get_node_attributes(graph, "atom_name").values())
        )
    return "CA" in atom_names


def extract_protein_sequence(graph):
    """Extract the protein sequence from a graph.

    This function extracts the protein sequence from the molecule represented
    by the input graph. It looks for the residue names corresponding to the C-alpha
    (CA) atoms in each node and converts them to single-letter amino acid codes.

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
    # graph.nodes.items() returns an iterable of tuples (node_id, node_attributes).
    # Examples:
    # (0, {'atom_name': 'N', 'residue_id': 1, 'residue_name': 'MET', 'label': 1})
    # (1, {'atom_name': 'CA', 'residue_id': 1, 'residue_name': 'MET', 'label': 5})
    # (2, {'atom_name': 'CB', 'residue_id': 1, 'residue_name': 'MET', 'label': 7})
    for _, node_attr in sorted(graph.nodes.items(), key=lambda x: x[0]):
        if node_attr["atom_name"] == "CA":
            protein_sequence.append(
                mol_def.AMINO_ACID_DICT.get(node_attr["residue_name"], "?")
            )
    info_seq["sequence"] = "".join(protein_sequence)
    info_seq["nb_res"] = len(protein_sequence)
    return info_seq


def export_protein_sequence_into_FASTA(protein_sequence_dict, filepath_name):
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
            logger.info(seq)
            content = f">Protein: {nb_res} residues\n" + "\n".join(seq)
            file.write(f"{content}\n")


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
    threshold = max(mol_def.BOND_LENGTH.values())
    logger.success(f"Bond threshold: {threshold} Angstrom")

    # molecular_system = read_structure_file_remove_hydrogens(input_file_path)
    molecular_system = remove_hydrogene(input_file_path)
    molecular_system, count_ion_solvant = count_remove_ion_solvant(
        molecular_system, input_file_path
    )

    atom_pairs = get_atom_pairs2(molecular_system, threshold)

    graph_return = convert_atom_pairs_to_graph(atom_pairs, molecular_system)

    graph_with_node_attributes = add_attributes_to_nodes(graph_return, molecular_system)
    graph_list = get_graph_components(graph_with_node_attributes)

    if check_overlapping_residue:
        check_overlapping_residue_between_graphs(graph_list)

    graph_count_dict = count_molecule(graph_list)

    for index_graph, graph in enumerate(graph_count_dict.keys(), start=1):
        print_graph_fingerprint(graph, index_graph)

    graph_count_dict.update(count_ion_solvant)
    print_graph_inventory(graph_count_dict)

    filename = Path(input_file_path).stem
    if draw_graph_option:
        logger.info("Drawing graphs...")
        filename = Path(input_file_path).stem
        for index_graph, graph_count in enumerate(graph_count_dict.keys()):
            if isinstance(graph_count, nx.Graph):
                print_graph(graph_count, f"{filename}_{index_graph}.png")

    protein_sequence_dict = {}
    for index_graph, graph in enumerate(graph_count_dict.keys(), start=1):
        if isinstance(graph, nx.Graph) and is_protein(graph):
            protein_sequence_dict[index_graph] = extract_protein_sequence(graph)

    export_protein_sequence_into_FASTA(protein_sequence_dict, f"{filename}.fasta")
    fsiPDB.main(f"{filename}.fasta")


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

    if filename.suffix not in (".gro", ".pdb"):
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
        "--drawgraph",
        help="Draw graph of each molecule. Default: False.",
        default=False,
        action="store_true",
    )
    parser.add_argument(
        "--checkoverlapping",
        help="Check if some residues are overlapping between residues. Default: False.",
        default=False,
        action="store_true",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_arg()
    main(
        args.input,
        draw_graph_option=args.drawgraph,
        check_overlapping_residue=args.checkoverlapping,
    )
