from __future__ import annotations
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from ._typing import PathLike, Json


MOLECULE_DEFINITIONS_DATABASE_PATH = (
    Path(__file__).parent / "data" / "databases" / "molecule_definitions.json"
)
assert MOLECULE_DEFINITIONS_DATABASE_PATH.is_file()


@dataclass
class SmallMoleculeDefinition:
    """Base class for a small molecule as represented in grodecoder molecule definition database.

    Attributes
    ----------

    name: str
        Descriptive name
    residue_name: str
        Ion residue name as expected in topology
    atom_names: list[str]
        List of the atom composing the ion
    """

    name: str
    residue_name: str
    atom_names: list[str]

    def __hash__(self):
        return hash((self.name, self.residue_name, tuple(self.atom_names)))

    def almost_eq(self, other) -> bool:
        """Returns True if `residue_name` and `atom_names` members are equal."""
        if not isinstance(other, SmallMoleculeDefinition):
            return False
        return self.residue_name == other.residue_name and set(self.atom_names) == set(other.atom_names)


class IonDefinition(SmallMoleculeDefinition):
    """Ion as represented in grodecoder molecule definition database.

    Examples
    --------

    >>> po4 = IonDefinition(
    ...     name="pO4 ion - CHARMM",
    ...     residue_name="PO4",
    ...     atom_names=["P", "O1", "O2", "O3", "O4"]
    ... )
    """


class SolventDefinition(SmallMoleculeDefinition):
    """Solvent molecule as represented in grodecoder molecule definition database.

    Examples
    --------
    >>> octadecane = SolventDefinition(
    ...     name="octadecane - CG model with MARTINI",
    ...     residue_name='OD',
    ...     atom_names=["C", "C", "C", "C", "C"]
    ... )
    """


class DuplicateDefinitionError(Exception):
    """Raised when the same molecule is defined multiple times in the database."""

    def __init__(self, SmallMoleculeDefinition: SmallMoleculeDefinition):
        self.definition = SmallMoleculeDefinition

    def __str__(self):
        return f"Duplicate definition found for {self.definition.residue_name} with atom names {self.definition.atom_names}"


def read_molecule_database(path: PathLike) -> dict[str, Json]:
    """Reads the grodecoder molecule database.

    The database is a JSON file that contains information about
    different molecules, including their residue names.

    Returns
    -------
    dict[str, Any]
        A dictionary containing the molecule database.
        The keys are ("ions", "solvents", "dna", "rna", "amino_acids", "lipids").
        The values depend on the type of molecule.
    """
    with open(path, "r") as f:
        source_database = json.load(f)

    # Safety checks and reformat.
    molecule_database = source_database.copy()

    for molecule_type in ("ions", "solvents"):
        molecule_database[molecule_type] = {}

        for molecule in source_database[molecule_type]:
            residue_name = molecule["residue_name"]
            molecule_database[molecule_type].setdefault(residue_name, [])

            # Check if the molecule is already defined.
            current_molecule = SmallMoleculeDefinition(**molecule)
            for existing_molecule in molecule_database[molecule_type][residue_name]:
                if current_molecule.almost_eq(existing_molecule):
                    raise DuplicateDefinitionError(current_molecule)

            # Add the entry to the database.
            molecule_database[molecule_type][residue_name].append(
                current_molecule)

    return molecule_database


MOLECULE_DEFINITIONS = None
if MOLECULE_DEFINITIONS is None:
    MOLECULE_DEFINITIONS = read_molecule_database(
        MOLECULE_DEFINITIONS_DATABASE_PATH)


def get_ion_definitions() -> list[SmallMoleculeDefinition]:
    return MOLECULE_DEFINITIONS["ions"]


def get_solvent_definitions() -> list[SmallMoleculeDefinition]:
    return MOLECULE_DEFINITIONS["solvents"]
