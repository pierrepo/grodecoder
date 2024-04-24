from Bio import SeqIO
from pathlib import Path

import json
from loguru import logger
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
        # "return_type": "polymer_entity",
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
        
    for index_sequence, values in enumerate(dict_IdPDB_seq.values(), start=1):
        PDB_ID, sequence = values.values()
        logger.info(f"FASTA sequence {index_sequence} (first 20AA): {sequence[:20]}")
        logger.info(f"Putative PDB IDs(first 10 match): {PDB_ID}")
