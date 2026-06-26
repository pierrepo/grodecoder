import json
import re
from pathlib import Path

import click


def parse_rtp(filepath):
    """
    Parse a single GROMACS .rtp file and extract residues with their atoms.

    Parameters
    ----------
    filepath : str or Path
        Path to the .rtp file to parse.

    Returns
    -------
    dict
        Dictionary mapping residue names to lists of atom names.
    """
    residues = {}
    current_residue = None
    in_atoms = False

    with open(filepath) as f:
        for line in f:
            line = line.strip()

            if not line or line.startswith(";"):
                continue

            match_residue = re.match(r"^\[\ *(\w+)\ *\]$", line)
            if match_residue:
                name = match_residue.group(1)

                if name not in ("atoms", "bonds", "angles", "dihedrals", "impropers", "exclusions"):
                    current_residue = name
                    residues[current_residue] = []
                    in_atoms = False
                elif name == "atoms" and current_residue:
                    in_atoms = True
                else:
                    in_atoms = False

                continue

            if in_atoms and current_residue:
                parts = line.split()
                if parts:
                    residues[current_residue].append(parts[0])

    return residues


def parse_all_ff(top_directory):
    """
    Parse all .rtp files in every force field directory under a top directory.

    Parameters
    ----------
    top_directory : str or Path
        Path to the GROMACS top/ directory containing .ff subdirectories.

    Returns
    -------
    dict
        Dictionary mapping residue names to a list of unique atom compositions.
    """
    result = {}
    top_path = Path(top_directory)

    for entry in top_path.iterdir():
        if not entry.is_dir() or not entry.name.endswith(".ff"):
            continue

        ff_name = entry.name
        print(f"Working on force field: {ff_name}")

        rtp_files_found = False

        for file in entry.iterdir():
            if not file.is_file() or not file.name.endswith(".rtp"):
                continue

            rtp_files_found = True
            print(f"  Reading: {file.name}")

            residues = parse_rtp(file)

            for residue, atoms in residues.items():
                if residue not in result:
                    result[residue] = []
                if atoms not in result[residue]:
                    result[residue].append(atoms)

        if not rtp_files_found:
            print(f"No .rtp files found in {entry}")

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
    """Parse GROMACS .rtp files and extract residues and their atoms."""
    print(f"Starting parsing of force fields in: {top_directory}")

    result = parse_all_ff(top_directory)

    print(f"{len(result)} unique residues found across all force fields")

    output_path = Path(output)
    with open(output_path, "w") as f:
        json.dump(result, f, indent=2)

    print(f"Results saved to {output_path}")


if __name__ == "__main__":
    main()
