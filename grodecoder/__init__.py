import MDAnalysis as mda

from . import main
from .databases import read_molecule_database
from .typing import PathLike, UniverseLike


__all__ = ["main", "read_molecule_database", "read_topology", "number_of_atoms", "number_of_residues"]



def read_topology(path: PathLike) -> mda.Universe:
    """Reads a topology file."""
    return mda.Universe(path)


def number_of_atoms(universe: UniverseLike) -> int:
    """Returns the number of atoms in the universe."""
    return len(universe.atoms)


def number_of_residues(universe: UniverseLike) -> int:
    """Returns the number of residues in the universe."""
    return len(universe.residues)
