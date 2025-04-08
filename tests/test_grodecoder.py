import json
import sys
from pathlib import Path

import pytest

# Adds project root directory to the system path.
# Necessary as long as grodecoder is not a proper package.
# TODO: make grodecoder a package.
parent_dir = Path(__file__).resolve().parents[1]
sys.path.append(str(parent_dir))

import grodecoder as gd


# Path to directory containing files required for tests.
TEST_DATA_ROOT = Path(__file__).resolve().parent / "data"

# Path to the package examples directory.
EXAMPLE_DIR = Path(__file__).resolve().parent.parent / "data/examples"


def test_path_pair(topology: str, result: str) -> (Path, Path):
    return EXAMPLE_DIR / topology, TEST_DATA_ROOT / result


# @pytest.mark.parametrize(
#     "path_structure_file, path_reference_json",
#     [
#
#         test_path_pair("1BRS.gro", "1BRS.json"),
#         test_path_pair("1QJ8_ETH_ACN_MET_URE_SOL.pdb", "1QJ8_PACKMOL.json"),
#         test_path_pair("1QJ8_membrane.gro", "1QJ8_membrane.json"),
#         test_path_pair("1QJ8_solution.gro", "1QJ8_solution.json"),
#         test_path_pair("4MQJ_ABCD.gro", "4MQJ_ABCD.json"),
#         test_path_pair("4ZRY.gro", "4ZRY.json"),
#         test_path_pair("5MBA.gro", "5MBA.json"),
#         test_path_pair("5ZOA.gro", "5ZOA.json"),
#         test_path_pair("barstar.gro", "barstar.json"),
#         test_path_pair("DMPC_PI.gro", "DMPC_PI.json"),
#         test_path_pair("noriega_AA_CRD_3CAL.gro", "data/3CAL_AA.json"),
#         test_path_pair("noriega_CG_CRD_3CAL.gro", "data/3CAL_CG.json"),
#         test_path_pair("2MAT.gro", "data/2MAT.json"),
#     ],
# )
# def test_with_param(path_structure_file, path_reference_json):
#     path_test_json = gd.main(path_structure_file)
#     with open(path_test_json, "r") as f:
#         JSON_content_test = json.load(f)
#
#     with open(path_reference_json, "r") as f:
#         JSON_content_data = json.load(f)
#
#     key_to_check = [
#         "id",
#         "number_of_atoms",
#         "number_of_molecules",
#         "residue_names",
#         "residue_ids",
#         "formula_without_h",
#         "molecular_type",
#         "sequence",
#         "comment",
#     ]
#
#     for index in range(len(JSON_content_data["inventory"])):
#         for key in key_to_check:
#             assert JSON_content_test["inventory"][index][key] == JSON_content_data["inventory"][index][key]
#
#     assert JSON_content_test["resolution"] == JSON_content_data["resolution"]
#     assert JSON_content_test["file_md5sum"] == JSON_content_data["file_md5sum"]
