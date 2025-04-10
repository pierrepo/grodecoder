import csv
import json
from typing import Any

from icecream import ic

from grodecoder import mol_def as mol_def

Entry = dict[str, Any]


def mol_def_to_entries(definitions: list[dict]) -> list[Entry]:
    entries = []
    for item in definitions:
        entries.append(
            {
                "name": item["name"],
                "residue_name": item["res_name"],
                "atom_names": item["atom_names"],
            }
        )
    return entries


def mol_def_nucleic_to_entries(nucleic: dict[str | str]) -> list[Entry]:
    entries = []
    names = {
        "A": "adenine",
        "C": "cytosine",
        "G": "guanine",
        "U": "uracil",
        "T": "thymine",
    }
    for key, value in nucleic.items():
        entries.append(
            {
                "name": names[value],
                "residue_name": key,
                "short_name": value,
            }
        )
    return entries


def amino_acids_to_json():
    amino_acids = mol_def.AMINO_ACID_DICT
    entries = []
    names = {
        "A": "alanine",
        "C": "cysteine",
        "D": "aspartic acid",
        "E": "glutamic acid",
        "F": "phenylalanine",
        "G": "glycine",
        "H": "histidine",
        "I": "isoleucine",
        "K": "lysine",
        "L": "leucine",
        "M": "methionine",
        "N": "asparagine",
        "P": "proline",
        "Q": "glutamine",
        "R": "arginine",
        "S": "serine",
        "T": "threonine",
        "V": "valine",
        "W": "tryptophan",
        "Y": "tyrosine",
    }
    for long_name, short_name in amino_acids.items():
        entries.append(
            {
                "name": names[short_name],
                "residue_name": long_name,
                "short_name": short_name,
            }
        )
    return entries


def lipid_MAD_to_json():
    csv_path = "grodecoder/data/databases/lipid_MAD.csv"
    entries = []
    with open(csv_path, "r") as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            entry = {
                "name": row["Name"],
                "alias": row["Alias"],
                "category": row["Category"],
                "link": row["Lien"],
            }
            entries.append(entry)
    return entries


def lipid_CHARMM_to_json():
    csv_path = "grodecoder/data/databases/lipid_CHARMM_GUI_CSML.csv"
    entries = []
    with open(csv_path, "r") as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            entry = {
                "name": row["Name"],
                "alias": row["Alias"],
                "category": row["Category"],
                "links": {
                    "view": row["View_Link"],
                    "pdb": row["PDB_Link"],
                },
                "formula": row["Formula"],
                "residue_name": row["Res_name_PDB"],
            }
            entries.append(entry)
    return entries


def ions_to_json():
    return mol_def_to_entries(mol_def.IONS_LIST)


def solvents_to_json():
    return mol_def_to_entries(mol_def.SOLVANTS_LIST)


def dna_to_json():
    return mol_def_nucleic_to_entries(mol_def.NUCLEIC_ACIDS_DNA)


def rna_to_json():
    return mol_def_nucleic_to_entries(mol_def.NUCLEIC_ACIDS_RNA)



def main():
    all_entries = {
        "ions": ions_to_json(),
        "solvents": solvents_to_json(),
        "dna": dna_to_json(),
        "rna": rna_to_json(),
        "amino_acids": amino_acids_to_json(),
        "lipids": {
            "MAD": lipid_MAD_to_json(),
            "CHARMM": lipid_CHARMM_to_json(),
        },
    }
    with open("molecule_definitions.json", "w") as f:
        json.dump(all_entries, f, indent=2)

if __name__ == "__main__":
    main()
