from Bio import SeqIO
from pathlib import Path

import json
import requests


def API_PDB_search_based_sequence(sequence, max_element=10):
    """This function searches the Protein Data Bank (PDB) database based on a given sequence and returns a list of PDB IDs.

    Ressources
    ----------
        https://education.molssi.org/python-scripting-biochemistry/chapters/rcsb_api.html
        https://search.rcsb.org/index.html#return-count
        
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
    my_query = {
        "query": {
            "type": "terminal",
            "service": "sequence",
            "parameters": {
                "evalue_cutoff": 0.1,
                "identity_cutoff": 0,
                "sequence_type": "protein",
                "value": sequence
            }
        },
        # "return_type": "polymer_entity", # Returns a list of PDB IDs appended with entity IDs in the format of a [pdb_id]_[entity_id], corresponding to polymeric molecular entities. 
        "return_type": "entry",
        "request_options": {
            "results_verbosity": "compact",
            "results_content_type": ["experimental"],
            "scoring_strategy": "combined",
            "sort": [{"sort_by": "score", "direction": "desc"}]
        }
    }

    my_query = json.dumps(my_query)
    data = requests.get(f"https://search.rcsb.org/rcsbsearch/v2/query?json={my_query}")
    results = data.json()
    return results["result_set"][:max_element]


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
