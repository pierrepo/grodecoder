from dataclasses import dataclass

import grodecoder as gd
import grodecoder.v1 as v1
import grodecoder.databases as DB

# DEBUG
import icecream

icecream.install()

# Import the "old" grodecoder for comparison.
# To be removed.
v1.logger.remove()


@dataclass
class SmallMoleculeCount:
    name: str
    residue_name: str
    number_of_atoms: int
    number_of_residues: int


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



# def count_ions_solvents_grodecoder(universe: UniverseLike) -> dict[str, int]:
#     warnings.filterwarnings("ignore")
#     _, counts = v1.count_remove_ion_solvant(universe, "foo")
#     result = {}
#     for molecule, dict_count in counts.items():
#         count = dict_count.get("graph")
#         res_name = " ".join(
#             set(nx.get_node_attributes(molecule, "residue_name").values())
#         )
#         result[res_name] = count
#     return result


# def create_inventory(ion_count: dict[str, SmallMoleculeCount], solvent_count: dict[str, SmallMoleculeCount]):
#     inventory = []
#     for count in ion_count.values():
#         item = {
#             "id": len(inventory) + 1,
#             "number_of_atoms": count.number_of_atoms,
#             "number_of_molecules": count.number_of_residues,
#             "residue_names": count.residue_name,
#         }
#         inventory.append(item)
#
#     for residue_name, count in solvent_count.items():
#         item = {
#             "id": len(inventory) + 1,
#             "number_of_atoms": count.number_of_atoms,
#             "number_of_molecules": count.number_of_residues,
#             "residue_names": count.residue_name,
#         }
#         inventory.append(item)


def main():
    topology_path = "grodecoder/data/examples/1BRS.gro"
    universe = gd.read_topology(topology_path)

    water = universe.select_atoms("resname TIP3")
    print(len(water))
    print(len(water.residues))

    ion_counts = count_solvents(universe)

    for i in ion_counts:
        ic(i)


if __name__ == "__main__":
    main()
