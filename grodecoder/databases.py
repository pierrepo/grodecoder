import json
from dataclasses import dataclass
from pathlib import Path

from ._typing import PathLike, Json


MOLECULE_DEFINITIONS_DATABASE_PATH = Path(__file__).parent / "data" / "databases" / "molecule_definitions.json"
assert MOLECULE_DEFINITIONS_DATABASE_PATH.is_file()


# TODO: convert to pydantic
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
        molecule_database = json.load(f)
    return molecule_database


MOLECULE_DEFINITIONS = None
if MOLECULE_DEFINITIONS is None:
    MOLECULE_DEFINITIONS = read_molecule_database(MOLECULE_DEFINITIONS_DATABASE_PATH)


def get_ion_definitions() -> list[IonDefinition]:
    # Convert to `IonDefinition`.
    # TODO: automate this using pydantic
    return [IonDefinition(**molecule) for molecule in MOLECULE_DEFINITIONS["ions"]]


def get_solvent_definitions() -> list[SolventDefinition]:
    # TODO: automate this using pydantic
    return [SolventDefinition(**molecule) for molecule in MOLECULE_DEFINITIONS["solvents"]]


