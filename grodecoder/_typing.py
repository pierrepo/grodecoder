"""Defines type aliases for the project."""

from pathlib import Path
import MDAnalysis as mda

PathLike = str | Path
UniverseLike = mda.core.universe.Universe | mda.core.groups.AtomGroup
Json = dict[str, "Json"] | list["Json"] | str | int | float | bool | None
