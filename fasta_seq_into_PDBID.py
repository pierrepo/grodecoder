from Bio import SeqIO
from pathlib import Path

# Reference: https://rcsbsearchapi.readthedocs.io/en/latest/quickstart.html#syntax
# Must execute "pip install rcsbsearchapi" before run this file
from rcsbsearchapi.search import SequenceQuery


def API_PDB_search_based_sequence(sequence, max_element=10):
    """This function searches the Protein Data Bank (PDB) database based on a given sequence and returns a list of PDB IDs.

    Parameters
    ----------
        sequence : str
            The amino acid sequence to search for in the PDB database.
        max_element : int, optional
            The maximum number of PDB IDs to retrieve. Defaults to 10.

    Returns
    -------
        list
            A list of PDB IDs corresponding to the sequences found in the PDB database.
    """
    list_ID_PDB = []
    # SequenceQuery(sequence, ..., sequence identity)
    results = SequenceQuery(sequence, 1, 0.9)
    
    # results("polymer_entity") produces an iterator of IDs with return type - polymer entities
    for index_result, polyid in enumerate(results("polymer_entity")):
        if index_result == max_element: break
        list_ID_PDB.append(polyid)
    return list_ID_PDB


def fasta_format_IdPDB(filepath, dict_IdPDB_seq):
    """This function formats the sequences corresponding to PDB IDs in a FASTA-like format and writes them to a file.

    Parameters
    ----------
        filepath : str
            The path to the output file.
        dict_IdPDB_seq : dict
            A dictionary containing PDB IDs as keys and their corresponding sequences as values.
    """
    with open(f"PDB_ID_{Path(filepath).stem}.{Path(filepath).suffix}", "w") as file:
        for values in dict_IdPDB_seq.values():
            IdPDB, sequence = values.values()
            
            # For only have 80 residues for each line
            seq = [sequence[i : i + 80] for i in range(0, len(sequence), 80)]
            content = f">First {len(IdPDB)} PDB ID: {' '.join(IdPDB)}\n" + "\n".join(seq)
            file.write(f"{content}\n")
        
            
def main(filepath):
    """ Based on sequence from the input FASTA file, search first 10 PDB id (by default) and put-it back into a FASTA file format.

    Parameters
    ----------
        filepath: str
            The path to the input FASTA file.
    """
    records = list(SeqIO.parse(filepath, "fasta"))
    sequences = [str(seq.seq) for seq in records]

    dict_IdPDB_seq = {}
    for index, sequence in enumerate(sequences):
        results = API_PDB_search_based_sequence(sequence)
        dict_IdPDB_seq[index] = {"IdPDB": results, 
                                "sequence": sequence}
    fasta_format_IdPDB(filepath, dict_IdPDB_seq)
