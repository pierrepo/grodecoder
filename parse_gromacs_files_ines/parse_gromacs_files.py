import json
import re
from pathlib import Path

import click
from loguru import logger
from pydantic import BaseModel


class ResidueEntry(BaseModel):
    name: str
    atoms: list[list[str]]
    molecular_type: str


class ResidueDatabase(BaseModel):
    residues: list[ResidueEntry]


RTP_FILENAME_TO_TYPE = {
    "aminoacids": "RESIDUE",
    "dna": "DNA",
    "rna": "RNA",
    "lipids": "LIPID",
}


def get_molecular_type_from_rtp_filename(filename: str) -> str:
    """
    get the molecular type from the name of a .rtp file.

    Parameters
    ----------
    filename : str
        Name of the .rtp file (e.g. 'aminoacids.rtp', 'dna.rtp').

    Returns
    -------
    str
        Molecular type string. Defaults to 'RESIDUE' if the filename is unknown.
    """
    stem = Path(filename).stem.lower()
    logger.info(f"stem : {stem}")
    return RTP_FILENAME_TO_TYPE.get(stem, "RESIDUE")


def get_molecular_type_from_itp_filename(filename: str) -> str | None:
    """
    Get the molecular type from the name of a .itp file.

    Parameters
    ----------
    filename : str
        Name of the .itp file (e.g. 'ions.itp', 'water.itp', 'ffnonbonded.itp').

    Returns
    -------
    str | None
        Molecular type string, or None if the file should be skipped.
    """
    stem = Path(filename).stem.lower()
    if stem.startswith("ff"):
        logger.info(f"stem : {stem}")
        return None
    if stem.startswith("ion"):
        logger.info(f"stem : {stem}")
        return "ION"
    return "WATER"


def parse_rtp(filepath: Path) -> dict[str, list[str]]:
    """
    Parse a single GROMACS .rtp file and extract residues with their atoms.

    Parameters
    ----------
    filepath : Path
        Path to the .rtp file to parse.

    Returns
    -------
    dict
        Mapping residue name -> list of atom names (first column of [ atoms ] section).
    """
    residues = {}
    current_residue = None
    in_atoms = False

    section_keywords = {"atoms", "bonds", "angles", "dihedrals", "impropers", "exclusions", "bondedtypes"}

    with open(filepath) as f:
        for line in f:
            line = line.strip()

            if not line or line.startswith(";"):
                continue

            match_section = re.match(r"^\[\ *(\w+)\ *\]$", line)
            if match_section:
                logger.success(f"match section : {match_section}")
                section_name = match_section.group(1)

                if section_name not in section_keywords:
                    current_residue = section_name
                    logger.info(f"Current residue : {current_residue}")
                    residues[current_residue] = []
                    in_atoms = False
                elif section_name == "atoms" and current_residue:
                    in_atoms = True
                else:
                    in_atoms = False
                continue

            if in_atoms and current_residue:
                parts = line.split()
                if parts:
                    logger.info(f"Added atom : {parts[0]}")
                    residues[current_residue].append(parts[0])

    return residues


def parse_itp(filepath: Path) -> dict[str, list[str]]:
    """
    Parse a single GROMACS .itp file and extract molecules with their atoms.

    Parameters
    ----------
    filepath : Path
        Path to the .itp file to parse.

    Returns
    -------
    dict
        Mapping molecule name -> list of atom names (column 5 of [ atoms ] section).
    """
    molecules = {}
    current_molecule = None
    in_moleculetype = False
    in_atoms = False

    with open(filepath) as f:
        for line in f:
            line = line.strip()

            if not line or line.startswith(";") or line.startswith("#"):
                continue

            match_section = re.match(r"^\[\ *(\w+)\ *\]$", line)
            if match_section:
                logger.success(f"match section : {match_section}")
                section_name = match_section.group(1)

                if section_name == "moleculetype":
                    in_moleculetype = True
                    in_atoms = False
                    current_molecule = None
                elif section_name == "atoms" and current_molecule:
                    in_atoms = True
                    in_moleculetype = False
                else:
                    in_atoms = False
                    in_moleculetype = False
                continue

            if in_moleculetype:
                parts = line.split()
                if parts:
                    current_molecule = parts[0]
                    logger.info(f"Current molecule: {current_molecule}")
                    molecules[current_molecule] = []
                    in_moleculetype = False
                continue

            if in_atoms and current_molecule:
                parts = line.split()
                if len(parts) >= 5:
                    molecules[current_molecule].append(parts[4])
                    logger.info(f"Added atom : {parts[4]}")

    return molecules


def parse_all_ff(top_directory: Path) -> ResidueDatabase:
    """
    Parse all .rtp and .itp files across every force field directory in top_directory.

    Parameters
    ----------
    top_directory : Path
        Path to the GROMACS top/ directory containing .ff subdirectories.

    Returns
    -------
    ResidueDatabase
        Pydantic model containing all parsed entries.
    """
    aggregated: dict[str, dict] = {}

    for entry in top_directory.iterdir():
        if not entry.is_dir() or not entry.name.endswith(".ff"):
            continue

        ff_name = entry.name
        logger.info(f"Working on force field: {ff_name}")

        for file in entry.iterdir():
            if not file.is_file() or file.suffix != ".rtp":
                continue

            logger.info(f"  Reading rtp: {file.name}")
            molecular_type = get_molecular_type_from_rtp_filename(file.name)
            residues = parse_rtp(file)

            for residue_name, atoms in residues.items():
                if residue_name not in aggregated:
                    aggregated[residue_name] = {
                        "compositions": [],
                        "molecular_type": molecular_type,
                    }
                atoms_tuple = tuple(atoms)
                if atoms_tuple not in aggregated[residue_name]["compositions"]:
                    aggregated[residue_name]["compositions"].append(atoms_tuple)

        for file in entry.iterdir():
            if not file.is_file() or file.suffix != ".itp":
                continue

            molecular_type = get_molecular_type_from_itp_filename(file.name)
            if molecular_type is None:
                logger.info(f"  Skipping itp: {file.name}")
                continue

            logger.info(f"  Reading itp: {file.name}")
            molecules = parse_itp(file)

            for molecule_name, atoms in molecules.items():
                if not atoms:
                    continue

                if molecule_name not in aggregated:
                    aggregated[molecule_name] = {
                        "compositions": [],
                        "molecular_type": molecular_type,
                    }
                atoms_tuple = tuple(atoms)
                if atoms_tuple not in aggregated[molecule_name]["compositions"]:
                    aggregated[molecule_name]["compositions"].append(atoms_tuple)

    if not aggregated:
        logger.info(f"No .ff directories found in {top_directory}")

    residue_entries = []
    for name, data in aggregated.items():
        residue_entries.append(
            ResidueEntry(
                name=name,
                atoms=[list(compo) for compo in data["compositions"]],
                molecular_type=data["molecular_type"],
            )
        )

    return ResidueDatabase(residues=residue_entries)


@click.command()
@click.argument(
    "top_directory",
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
)
@click.option(
    "--output",
    required=True,
    type=click.Path(),
    help="Path to the output JSON file.",
)
def main(top_directory, output):
    logger.info(f"Starting parsing of force fields in: {top_directory}")

    top_path = Path(top_directory)
    database = parse_all_ff(top_path)

    logger.info(f"{len(database.residues)} unique entries found across all force fields")

    output_path = Path(output)
    with open(output_path, "w") as f:
        json.dump(database.model_dump(), f, indent=2)

    logger.info(f"Results saved to {output_path}")


if __name__ == "__main__":
    main()
