import argparse
from Bio import SeqIO
import json
from pathlib import Path
import requests

from loguru import logger


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
    logger.info("Searching PDB IDs based on a given sequence")
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
        "return_type": "polymer_entity",
        # "return_type": "entry",
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


def get_macromolecular_names(PDB_ID, polymer_entity_id, set_names):
    """Retrieve macromolecular names associated with a given PDB ID and polymer entity ID.

    References
    ----------
    https://www.rcsb.org/docs/programmatic-access/web-services-overview#data-api
    https://data.rcsb.org/index.html#data-api
    https://data.rcsb.org/redoc/index.html#tag/Entity-Service/operation/getPolymerEntityById 

    Parameters
    ----------
        PDB_ID : str
            The PDB ID of the structure.
        polymer_entity_id : str
            The polymer entity ID within the PDB structure.
        set_names : set
            A set containing existing macromolecular names. Retrieved names will be added to this set.

    Returns
    -------
        set
            A set containing the updated list of macromolecular names.

    Raises
    ------
        requests.exceptions.RequestException
            If an error occurs while making the API request.
    """
    logger.info("Retrieving macromolecular names for one PDB ID")
    try : 
        my_query = requests.get(f"https://data.rcsb.org/rest/v1/core/polymer_entity/{PDB_ID}/{polymer_entity_id}")
        my_query.raise_for_status()
        
        results = my_query.json()
        rcsb_macromolecular_names_combined = results["rcsb_polymer_entity"]["rcsb_macromolecular_names_combined"]
        for index in range(len(rcsb_macromolecular_names_combined)):
            result = results["rcsb_polymer_entity"]["rcsb_macromolecular_names_combined"][index]["name"].upper()
            set_names.add(result)
        return set_names
    
    except requests.exceptions.RequestException as err:
        print(f"Error: {err}")
        return set_names
            

def treat_PDB_ID_to_macromolecular_names(PDB_ID_polymer_entity_id):
    """Retrieve macromolecular names for multiple PDB IDs and polymer entity IDs.

    Parameters
    ----------
        PDB_ID_polymer_entity_id : list of str
            A list containing strings in the format '(PDB_ID)_(polymer_entity_id)'.

    Returns
    -------
        set
            A set containing unique macromolecular names extracted from the specified PDB IDs and polymer entity IDs.
    """
    set_macromolecular_names = set()
    
    for ID_entity in PDB_ID_polymer_entity_id:
        PDB_ID, polymer_entity_id = ID_entity.split("_")
        get_macromolecular_names(PDB_ID, polymer_entity_id, set_macromolecular_names) 
    return set_macromolecular_names


def main(filepath):
    """ Based on sequence from the input FASTA file, search and print the first 10 PDB id (by default) and their macromolecular names.

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
        set_macromolecular_names = treat_PDB_ID_to_macromolecular_names(results)
        dict_IdPDB_seq[index] = {"PDB_ID": results, 
                                "sequence": sequence, 
                                "macromolecular_names": set_macromolecular_names}

    for index_sequence, values in enumerate(dict_IdPDB_seq.values(), start=1):
        PDB_ID, sequence, macromolecular_names = values.values()
        logger.success(f"FASTA sequence {index_sequence} (first 20 AA): {sequence[:20]}")
        logger.success(f"Putative macromolecular names: {macromolecular_names}")
        logger.success(f"Putative PDB IDs (first 10 match): {PDB_ID}")


def is_a_fasta_file(filepath):
    """Check if the given filepath points to an existing fasta file.

    Parameters
    ----------
        filepath : str
            Path of the file.

    Raises
    ------
        argparse.ArgumentTypeError
            If the given filepath is not an existing file,
            or if it does not have a fasta extension.

    Returns
    -------
        str
            The validated path.
    """
    filename = Path(filepath)
    if not Path.is_file(filename):
        raise argparse.ArgumentTypeError(f"{filepath} does not exist")

    if filename.suffix != ".fasta":
        raise argparse.ArgumentTypeError(f"{filepath} is not a .fasta file")
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
        prog="search_into_PDB",
        description="Seach PDB IDs and macromolecular names in PDB database, based on an extracted sequence from a fasta file.",
        usage="search_into_PDB.py [-h] --fasta fasta_file",
    )
    parser.add_argument(
        "--fasta",
        type=is_a_fasta_file,
        help="fasta file path",
        required=True,
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_arg()
    main(
        args.fasta
    )