import json
import os
from pathlib import Path
import sys

import pytest

parent_dir = Path(__file__).resolve().parents[1]
sys.path.append(str(parent_dir))

import grodecoder as gd


@pytest.mark.parametrize(
    "path_structure_file, path_reference_json",
    [
        (os.path.join(parent_dir, "data/examples/1BRS.gro"), "data/1BRS.json"),
        (
            os.path.join(parent_dir, "data/examples/1QJ8_ETH_ACN_MET_URE_SOL.pdb"),
            "data/1QJ8_PACKMOL.json",
        ),
        (
            os.path.join(parent_dir, "data/examples/1QJ8_membrane.gro"),
            "data/1QJ8_membrane.json",
        ),
        (
            os.path.join(parent_dir, "data/examples/1QJ8_solution.gro"),
            "data/1QJ8_solution.json",
        ),
        (
            os.path.join(parent_dir, "data/examples/4MQJ_ABCD.gro"),
            "data/4MQJ_ABCD.json",
        ),
        (os.path.join(parent_dir, "data/examples/4ZRY.gro"), "data/4ZRY.json"),
        (os.path.join(parent_dir, "data/examples/5MBA.gro"), "data/5MBA.json"),
        (os.path.join(parent_dir, "data/examples/5ZOA.gro"), "data/5ZOA.json"),
        (os.path.join(parent_dir, "data/examples/barstar.gro"), "data/barstar.json"),
        (os.path.join(parent_dir, "data/examples/DMPC_PI.gro"), "data/DMPC_PI.json"),
        (
            os.path.join(parent_dir, "data/examples/noriega_AA_CRD_3CAL.gro"),
            "data/3CAL_AA.json",
        ),
        (
            os.path.join(parent_dir, "data/examples/noriega_CG_CRD_3CAL.gro"),
            "data/3CAL_CG.json",
        ),
        (os.path.join(parent_dir, "data/examples/2MAT.gro"),"data/2MAT.json",),
    ],
)
def test_with_param(path_structure_file, path_reference_json):
    path_test_json = gd.main(path_structure_file)
    with open(path_test_json, "r") as f:
        JSON_content_test = json.load(f)

    with open(path_reference_json, "r") as f:
        JSON_content_data = json.load(f)

    key_to_check = [
        "id",
        "number_of_atoms",
        "number_of_molecules",
        "residue_names",
        "residue_ids",
        "formula_without_h",
        "molecular_type",
        "sequence",
        "comment",
    ]

    for index in range(len(JSON_content_data["inventory"])):
        for key in key_to_check:
            assert (
                JSON_content_test["inventory"][index][key]
                == JSON_content_data["inventory"][index][key]
            )

    assert JSON_content_test["resolution"] == JSON_content_data["resolution"]
    assert JSON_content_test["file_md5sum"] == JSON_content_data["file_md5sum"]
