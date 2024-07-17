"""GroDecoder extracts and identifies the molecular components of a structure file 
(PDB, GRO, .coor or .crd) issued from a molecular dynamics simulation.
And return a molecular inventory, in a JSON file, with information the composition: 
molecular type (protein, ion, solvant, nucleic acids or unknown), the number of atom, 
their occurence, their name, residues names, ...


Example
-------
Running the script after isntallation
    $ 

"""

__authors__ = "Karine Duong, Pierre Poulain"
__license__ = "BSD 3-Clause License"
__version__ = "1.1.0"
__maintainer__ = "Pierre Poulain"
__email__ = "pierre.poulain@cupnet.net"


from .grodecoder import (get_distance_matrix_between_atom, 
                         get_atom_pairs, 
                         get_atom_pairs_from_threshold, 
                         get_atom_pairs_from_guess_bonds, 
                         convert_atom_pairs_to_graph, 
                         add_attributes_to_nodes, 
                         get_graph_components, 
                         get_graph_fingerprint, 
                         get_graph_fingerprint_str, 
                         print_graph_fingerprint, 
                         get_intervals, 
                         extract_interval, 
                         get_formula_based_atom_name, 
                         count_molecule, 
                         print_graph_inventory, 
                         print_first_atoms, 
                         read_structure_file_remove_hydrogens, 
                         remove_hydrogene, 
                         find_ion_solvant, 
                         count_remove_ion_solvant, 
                         find_lipids, 
                         count_remove_lipid, 
                         check_overlapping_residue_between_graphs, 
                         is_protein, 
                         extract_protein_sequence, 
                         export_protein_sequence_into_FASTA, 
                         is_lipid, 
                         guess_resolution, 
                         get_git_last_commit_date, 
                         get_git_last_commit_hash, 
                         export_inventory, 
                         is_met, 
                         is_nucleic_acid, 
                         extract_nucleic_acid_sequence, 
                         main, 
                         is_a_structure_file, 
                         is_a_valid_threshold, 
                         parse_arg,
                         )


from .mol_def import (BOND_LENGTH, 
                      AMINO_ACID_DICT, 
                      NUCLEIC_ACIDS_DNA, 
                      NUCLEIC_ACIDS_RNA, 
                      NUCLEIC_ACIDS, 
                      IONS_LIST, 
                      SOLVANTS_LIST,
                      )


from .search_into_PDB import (API_PDB_search_based_sequence, 
                              get_macromolecular_names, 
                              treat_PDB_ID_to_macromolecular_names, 
                              get_macromolecular_names_bis, 
                              get_organism_names, 
                              get_info_one_pdb_id, 
                              main, 
                              is_a_fasta_file, 
                              parse_arg, 
                              )


__all__ = ["get_distance_matrix_between_atom", 
            "get_atom_pairs", 
            "get_atom_pairs_from_threshold", 
            "get_atom_pairs_from_guess_bonds", 
            "convert_atom_pairs_to_graph", 
            "add_attributes_to_nodes", 
            "get_graph_components", 
            "get_graph_fingerprint", 
            "get_graph_fingerprint_str", 
            "print_graph_fingerprint", 
            "get_intervals", 
            "extract_interval", 
            "get_formula_based_atom_name", 
            "count_molecule", 
            "print_graph_inventory", 
            "print_first_atoms", 
            "read_structure_file_remove_hydrogens", 
            "remove_hydrogene", 
            "find_ion_solvant", 
            "count_remove_ion_solvant", 
            "find_lipids", 
            "count_remove_lipid", 
            "check_overlapping_residue_between_graphs", 
            "is_protein", 
            "extract_protein_sequence", 
            "export_protein_sequence_into_FASTA", 
            "is_lipid", 
            "guess_resolution", 
            "get_git_last_commit_date", 
            "get_git_last_commit_hash", 
            "export_inventory", 
            "is_met", 
            "is_nucleic_acid", 
            "extract_nucleic_acid_sequence", 
            "main", 
            "is_a_structure_file", 
            "is_a_valid_threshold", 
            "parse_arg",
            "BOND_LENGTH", 
            "AMINO_ACID_DICT", 
            "NUCLEIC_ACIDS_DNA", 
            "NUCLEIC_ACIDS_RNA", 
            "NUCLEIC_ACIDS", 
            "IONS_LIST", 
            "SOLVANTS_LIST",
            "API_PDB_search_based_sequence", 
            "get_macromolecular_names", 
            "treat_PDB_ID_to_macromolecular_names", 
            "get_macromolecular_names_bis", 
            "get_organism_names", 
            "get_info_one_pdb_id", 
            "is_a_fasta_file",
           ]