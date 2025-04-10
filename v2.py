import warnings
from dataclasses import dataclass

import MDAnalysis as mda
import networkx as nx
from icecream import ic

from grodecoder import main as grodecoder
from grodecoder import mol_def as mol_def


Universe = mda.AtomGroup | mda.Universe


grodecoder.logger.remove()

@dataclass
class SmallMoleculeCount:
    residue_name: str
    number_of_atoms: int
    number_of_residues: int


def count(universe: Universe, definitions: list[dict]) -> dict[str, int]:
    count = {}
    for molecule_definition in definitions:
        residue_name = molecule_definition["res_name"]
        atom_names = molecule_definition["atom_names"]

        selection_str = f"resname {residue_name} and name {' '.join(atom_names)}" 
        selection = universe.select_atoms(selection_str)
        if len(selection) > 0:
            count[residue_name] = SmallMoleculeCount(
                residue_name=residue_name,
                number_of_atoms=len(selection),
                number_of_residues=len(selection.residues),
            )
    return count


def count_ions(universe: Universe) -> dict[str, int]:
    definitions = mol_def.IONS_LIST
    return count(universe, definitions)


def count_solvents(universe: Universe) -> dict[str, int]:
    definitions = mol_def.SOLVANTS_LIST
    return count(universe, definitions)


def count_ions_solvents(universe: Universe) -> dict[str, int]:
    definitions = mol_def.IONS_LIST + mol_def.SOLVANTS_LIST
    return count(universe, definitions)


def count_ions_solvents_grodecoder(universe: Universe) -> dict[str, int]:
    warnings.filterwarnings("ignore")
    _, counts = grodecoder.count_remove_ion_solvant(universe, "foo")
    result = {}
    for molecule, dict_count in counts.items():
        count = dict_count.get("graph")
        res_name = " ".join(
            set(nx.get_node_attributes(molecule, "residue_name").values())
        )
        result[res_name] = count
    return result


def create_inventory(ion_count: dict[str, SmallMoleculeCount], solvent_count: dict[str, SmallMoleculeCount]):
    inventory = []
    for count in ion_count.values():
        item = {
            "id": len(inventory) + 1,
            "number_of_atoms": count.number_of_atoms,
            "number_of_molecules": count.number_of_residues,
            "residue_names": count.residue_name,
        }
        inventory.append(item)

    for residue_name, count in solvent_count.items():
        item = {
            "id": len(inventory) + 1,
            "number_of_atoms": count.number_of_atoms,
            "number_of_molecules": count.number_of_residues,
            "residue_names": count.residue_name,
        }
        inventory.append(item)


def main():
    # topology = "grodecoder/data/examples/barstar.gro"
    topology = "grodecoder/data/examples/1BRS.gro"
    universe = mda.Universe(topology)

    ions = count_ions(universe)
    solvents = count_solvents(universe)


    exit()



if __name__ == "__main__":
    main()
