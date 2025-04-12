from dataclasses import dataclass
from loguru import logger

import grodecoder as gd
import grodecoder.v1 as v1
import grodecoder.databases as DB

# DEBUG
import icecream
import time

icecream.install()


# Import the "old" grodecoder for comparison.
# To be removed.
v1.logger.remove()


@dataclass
class SmallMoleculeCount:
    name: str
    residue_name: str
    number_of_atoms: int = 0
    number_of_residues: int = 0


def count_small_molecule(
    universe: gd.UniverseLike, definition: DB.SmallMoleculeDefinition
) -> SmallMoleculeCount:
    """Counts the number of a given molecule in a ensemble of atoms."""
    selection_str = f"resname {definition.residue_name} and name {' '.join(definition.atom_names)}"
    selection = universe.select_atoms(selection_str)

    return SmallMoleculeCount(
        name=definition.name,
        residue_name=definition.residue_name,
        number_of_atoms=gd.number_of_atoms(selection),
        number_of_residues=gd.number_of_residues(selection),
    )


def count_ions(universe: gd.UniverseLike) -> list[SmallMoleculeCount]:
    ions = DB.get_ion_definitions()
    counts = [count_small_molecule(universe, ion) for ion in ions]
    return [count for count in counts if count.number_of_atoms > 0]


def count_solvents(universe: gd.UniverseLike) -> list[SmallMoleculeCount]:
    solvents = DB.get_solvent_definitions()
    counts = [count_small_molecule(universe, solvent) for solvent in solvents]
    return [count for count in counts if count.number_of_atoms > 0]


class NoDefinitionFoundError(Exception):
    """Base class for all definition not found errors."""

    def __init__(self, residue_name: str, number_of_atoms: int):
        self.residue_name = residue_name
        self.number_of_atoms = number_of_atoms


class NoDefinitionWithMatchingNumberOfAtoms(NoDefinitionFoundError):
    """Exception raised when no definition is found for a given residue name."""

    def __str__(self):
        return f"No definition found for residue '{self.residue_name}' with {self.number_of_atoms} atoms"


class NoDefinitionWithMatchingAtomNames(NoDefinitionFoundError):
    """Exception raised when no definition is found for a given residue name."""

    def __str__(self):
        return f"No definition found for residue '{self.residue_name}' with atom names {self.atom_names}"


class MultipleDefinitionsWithMatchingAtomNames(NoDefinitionFoundError):
    """Exception raised when multiple definitions are found for a given residue name."""

    def __str__(self):
        return (
            f"Multiple definitions found for residue '{self.residue_name}' with atom names {self.atom_names}"
        )


def count2(universe) -> SmallMoleculeCount:
    solvents = DB.get_solvent_definitions()
    counts = {}

    for residue_name, definitions in solvents.items():
        selection = universe.select_atoms(f"resname {residue_name}")

        for residue in selection.residues:
            # Definitions with matching residue name and number of atoms.
            matches = [d for d in definitions if len(residue.atoms) == len(d.atom_names)]
            if not matches:
                raise NoDefinitionWithMatchingNumberOfAtoms(residue_name, len(residue.atoms))

            # Look for a definition with the same atom names.
            match = None
            for definition in matches:
                if set(definition.atom_names) == set(residue.atoms.names):
                    match = definition
                    break

            if match is None:
                raise NoDefinitionWithMatchingAtomNames(residue_name, residue.atoms.names)

            counts.setdefault(match, SmallMoleculeCount(definition.name, residue_name))
            counts[match].number_of_atoms += len(residue.atoms)
            counts[match].number_of_residues += 1
    return list(counts.values())


def main():
    topology_path = "grodecoder/data/examples/1BRS.gro"
    universe = gd.read_topology(topology_path)

    water = universe.select_atoms("resname TIP3")
    expected_number_of_atoms = len(water)
    expected_number_of_residues = len(water.residues)

    start = time.perf_counter()
    result = count2(universe)
    print(f"method 2: {time.perf_counter() - start:.2f} seconds")
    ic(result)


if __name__ == "__main__":
    main()
