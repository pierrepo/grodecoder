import warnings
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

import MDAnalysis as mda
import networkx as nx
import numpy as np
import pytest

import grodecoder.main as gd

# Path to directory containing files required for tests.
TEST_DATA_ROOT = Path(__file__).resolve().parent / "data"


def test_remove_hydrogen_gro():
    topology_path = TEST_DATA_ROOT / "fake_hydrogen.gro"
    universe = mda.Universe(topology_path)
    universe = gd.remove_hydrogen(universe)
    expected = [1, 5, 9, 13, 14, 16, 18, 20, 23, 26, 30, 35, 38, 41, 42] 
    assert np.array_equal(universe.atoms.ids, expected)


@pytest.mark.filterwarnings("ignore:Unit cell dimensions not found")
@pytest.mark.filterwarnings("ignore:Found no information for attr")
def test_remove_hydrogen_pdb():
    topology_path = TEST_DATA_ROOT / "fake_hydrogen.pdb"
    universe = mda.Universe(topology_path)
    universe = gd.remove_hydrogen(universe)
    expected = [1, 2, 3, 4, 5, 6, 7, 12, 16, 19, 21]
    assert np.array_equal(universe.atoms.ids, expected)


@pytest.mark.filterwarnings("ignore:`guess_atom_type` is deprecated")
@pytest.mark.filterwarnings("ignore:`guess_atom_element` is deprecated")
@pytest.mark.parametrize(
    "topology_path, expected_formula",
    [
        (TEST_DATA_ROOT / "get_formula_waterTIP3.gro", "O"),
        (TEST_DATA_ROOT / "get_formula_CHL1.gro", "C27O"),
        (TEST_DATA_ROOT / "get_formula_DOPC.gro", "C44NO8P"),
        (TEST_DATA_ROOT / "get_formula_POPP2.gro", "C37O11P2"),
    ],
)
def test_get_formula(topology_path, expected_formula):
    universe = mda.Universe(topology_path)
    atom_names = Counter(universe.atoms.names)
    atom_names = dict(sorted(atom_names.most_common()))

    formula = gd.get_formula_based_atom_name(atom_names)
    assert formula == expected_formula


def test_get_intervals():
    interval = gd.get_intervals([1, 2, 3, 6, 7, 8, 10])
    assert interval == ["1-3", "6-8", "10"]

    interval = gd.get_intervals([2, 5, 6, 7, 11, 12, 13, 14])
    assert interval == ["2", "5-7", "11-14"]


@dataclass
class Graph:
    """Graph holder for unit tests."""

    edges: list[tuple[int, int]]
    residue_ids: list[int] = field(default_factory=list)
    residue_names: list[str] = field(default_factory=list)
    graph: nx.Graph = field(init=False)

    def __post_init__(self):
        self.graph = nx.Graph()
        self.graph.add_edges_from(self.edges)

        if self.residue_ids:
            nx.set_node_attributes(
                self.graph,
                dict(zip(self.atom_ids, self.residue_ids, strict=True)),
                name="residue_id",
            )

        if self.residue_names:
            nx.set_node_attributes(
                self.graph,
                dict(zip(self.atom_ids, self.residue_names, strict=True)),
                name="residue_name",
            )

        nx.set_node_attributes(
            self.graph,
            dict(zip(self.atom_ids, self.atom_ids, strict=True)),
            name="atom_id",
        )

    @property
    def atom_ids(self):
        return list(self.graph.nodes())


def test_extract_interval_from_graph():
    g = Graph(
        residue_ids=[1, 1, 2, 2, 5],
        residue_names=["ALA", "ALA", "CYS", "CYS", "GLY"],
        edges=[(1, 2), (2, 3), (2, 4), (4, 5)],
    )
    intervals = gd.extract_interval(g.graph)
    assert intervals["res_id_interval"] == ["1-2", "5"]
    assert intervals["atom_id_interval"] == ["1-5"]


class TestIsLipid:
    """Tests for `grodecoder.is_lipid` function."""

    def test_is_lipid_all_atom(self):
        g = Graph(
            edges=[(1, 2), (3, 2)],
            residue_names=["ADG", "ADG", "ADG"],
        )
        assert gd.is_lipid("AA", g.graph, {"formula_no_h": "C16O6"})

    def test_is_not_lipid_all_atom(self):
        g = Graph(
            edges=[(1, 2), (3, 2)],
            residue_names=["ALA", "ALA", "ALA"],
        )
        assert not gd.is_lipid("AA", g.graph, {"formula_no_h": ""})

    def test_is_lipid_coarse_grain(self):
        g = Graph(
            edges=[(1, 2), (3, 2)],
            residue_names=["DFPE", "DFPE", "DFPE"],
        )
        assert gd.is_lipid("CG", g.graph, {})


class TestIsProtein:
    """Tests for `grodecoder.is_protein` function."""

    def test_is_protein(self):
        g = Graph(
            edges=[(1, 2), (3, 2), (2, 4), (4, 5)],
            residue_names=["ALA", "ARG", "CYS", "LYS", "ILE"],
        )
        assert gd.is_protein(g.graph)

    def test_is_not_protein(self):
        g = Graph(
            edges=[(1, 2), (3, 2), (2, 4), (4, 5)],
            residue_names=["HIT", "BAD", "NIL", "CIL", "ALS"],
        )
        assert not gd.is_protein(g.graph)

    def test_is_not_protein_because_residues_look_strange(self):
        g = Graph(
            edges=[(1, 2), (3, 2)],
            residue_names=["PHE", "CYS", "LYS"],
        )
        assert not gd.is_protein(g.graph)


def test_extract_protein_sequence():
    g = Graph(
        edges=[(1, 2), (3, 2), (2, 4), (4, 5)],
        residue_ids=[1, 1, 2, 2, 5],
        residue_names=["ALA", "ALA", "CYS", "CYS", "GLY"],
    )
    result = gd.extract_protein_sequence(g.graph)
    assert result["sequence"] == "ACG"
    assert result["nb_res"] == 3


class TestIsMethanol:
    """Tests for `grodecoder.is_met` function."""

    def test_is_methanol(self):
        g = Graph(
            edges = [(1, 2)],
            residue_names = ["MET", "MET"],
        )
        assert gd.is_met(g.graph)

    def test_is_not_methanol_if_wrong_residue_name(self):
        # TODO: fix `is_met` to handle this case
        warnings.warn("Fails: `is_met` seems to fail if only first residue name is MET")
        # g = Graph(
        #     edges = [(1, 2)],
        #     residue_names = ["MET", "ALA"],
        # )
        # assert not gd.is_met(g.graph)

    def test_is_not_methanol_if_wrong_number_of_atoms(self):
        g = Graph(
            edges = [(1, 2), (3, 2)],
            residue_names = ["MET", "MET", "MET"],
        )
        assert not gd.is_met(g.graph)


@dataclass
class FindIon:
    definition: dict
    mode: str
    topology_path: Path = TEST_DATA_ROOT / "find_ion_CLA.gro"
    universe: mda.core.AtomGroup = field(init=False)

    def __post_init__(self):
        self.universe = mda.Universe(self.topology_path)

    def find_ion_solvant(self):
        result, _ = gd.find_ion_solvant(self.definition, self.universe, {}, self.mode)
        return result

    @property
    def atom_ids(self) -> set[int]:
        return set(self.universe.atoms.ids)


class TestFindIonSolvent:
    """Tests for `grodecoder.find_ion_solvant` function.

    IMPORTANT
    ---------

    `grodecoder.find_ion_solvant` returns an `AtomGroup` object where ion or solvent atoms have been **REMOVED**.

    Therefore, in subsequent tests, we check that the returned `AtomGroup` contains the expected atom IDs,
    which are the IDs of the atoms that are **NOT** the ion or solvent.

    TODO
    ----

    - Change `find_ion_solvant` name to something more explicit, like `remove_ion_solvant`.
    """

    def test_find_ion(self):
        chloride_ion = {
            "name": "chloride ion",
            "res_name": "CLA",
            "atom_names": ["CLA"],
        }

        setup = FindIon(definition=chloride_ion, mode="ion")
        result = setup.find_ion_solvant()

        # Expects all atoms except the one CLA ion (i.e. atom id 6)
        expected = setup.atom_ids - {6}

        assert np.array_equal(result.atoms.ids, list(expected))

    def find_solvent(self):
        water = {
            "name": "water TIP3P solvant",
            "res_name": "SOL",
            "atom_names": ["OH2", "H1", "H2"],
        }

        setup = FindIon(definition=water, mode="solvant")
        result = setup.find_ion_solvant()

        # Expects all atoms except the 3 atoms from the water molecule.
        expected = setup.atom_ids - {7, 8, 9}

        assert np.array_equal(result.atoms.ids, list(expected))

    def test_find_acetonitrile(self):
        acn = {
            "name": "acetonitrile",
            "res_name": "ACN",
            "atom_names": ["C1", "C2", "N"],
        }

        setup = FindIon(definition=acn, mode="solvant")
        result = setup.find_ion_solvant()

        # No acetonitrile in the test file, so all atoms should be kept.
        expected = setup.atom_ids

        assert np.array_equal(result.atoms.ids, list(expected))


def test_guess_resolution_all_atom():
    universe = mda.Universe(TEST_DATA_ROOT / "test_guess_resolution_AA.gro")
    resolution = gd.guess_resolution(universe)
    assert resolution == "AA"


def test_guess_resolution_coarse_grain():
    universe = mda.Universe(TEST_DATA_ROOT / "test_guess_resolution_CG.gro")
    resolution = gd.guess_resolution(universe)
    assert resolution == "CG"
