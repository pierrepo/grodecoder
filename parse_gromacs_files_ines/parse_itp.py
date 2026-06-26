import json
import re
from pathlib import Path

import click


def parse_itp(filepath):
    """
    Parse a single GROMACS .itp file and extract molecules with their atoms.

    Parameters
    ----------
    filepath : str or Path
        Path to the .itp file to parse.

    Returns
    -------
    dict
        Dictionary mapping molecule names to lists of atom names.
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
                name = match_section.group(1)
                if name == "moleculetype":
                    in_moleculetype = True
                    in_atoms = False
                    current_molecule = None
                elif name == "atoms" and current_molecule:
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
                    molecules[current_molecule] = []
                    in_moleculetype = False
                continue

            if in_atoms and current_molecule:
                parts = line.split()
                if len(parts) >= 5:
                    molecules[current_molecule].append(parts[4])

    return molecules


def parse_all_ff_itp(top_directory):
    """
    Parse all .itp files in every force field directory under a top directory.

    Parameters
    ----------
    top_directory : str or Path
        Path to the GROMACS top/ directory containing .ff subdirectories.

    Returns
    -------
    dict
        Dictionary mapping molecule names to a list of unique atom compositions.
    """
    result = {}
    top_path = Path(top_directory)

    for entry in top_path.iterdir():
        if not entry.is_dir() or not entry.name.endswith(".ff"):
            continue

        print(f"Working on force field: {entry.name}")

        itp_files_found = False
        for file in entry.iterdir():
            if not file.is_file() or not file.name.endswith(".itp"):
                continue

            itp_files_found = True
            print(f"  Reading: {file.name}")

            molecules = parse_itp(file)

            for molecule, atoms in molecules.items():
                if molecule not in result:
                    result[molecule] = []
                if atoms not in result[molecule]:
                    result[molecule].append(atoms)

        if not itp_files_found:
            print(f"  No .itp files found in {entry}")

    if not result:
        print(f"No .ff directories found in {top_directory}")

    return result


@click.command()
@click.option(
    "--top-directory",
    required=True,
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    help="Path to the GROMACS top/ directory containing .ff subdirectories.",
)
@click.option(
    "--output",
    required=True,
    type=click.Path(),
    help="Path to the output JSON file.",
)
def main(top_directory, output):
    """Parse GROMACS .itp files and extract molecules and their atoms."""
    print(f"Starting parsing of force fields in: {top_directory}")

    result = parse_all_ff_itp(top_directory)

    print(f"{len(result)} unique molecules found across all force fields")

    output_path = Path(output)
    with open(output_path, "w") as f:
        json.dump(result, f)

    print(f"Results saved to {output_path}")


if __name__ == "__main__":
    main()
