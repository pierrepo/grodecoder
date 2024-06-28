import json
import os
import sys


parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
print(parent_dir)
sys.path.append(parent_dir)

import grodecoder as gd
import pytest


@pytest.mark.parametrize(
        "filepath_test, JSON_filepath_data", 
        [ (os.path.join(parent_dir, "data/examples/1BRS.gro"), "data/1BRS.json"), 
          (os.path.join(parent_dir, "data/examples/1QJ8_ETH_ACN_MET_URE.pdb"), "data/1QJ8_PACKMOL.json"), 
          (os.path.join(parent_dir,"data/examples/1QJ8_membrane.gro"), "data/1QJ8_membrane.json"), 
          (os.path.join(parent_dir,"data/examples/1QJ8_solution.gro"), "data/1QJ8_solution.json"), 
          (os.path.join(parent_dir,"data/examples/4MQJ_ABCD.gro"), "data/4MQJ_ABCD.json"), 
          (os.path.join(parent_dir,"data/examples/4ZRY.gro"), "data/4ZRY.json"), 
          (os.path.join(parent_dir,"data/examples/5MBA.gro"), "data/5MBA.json"), 
          (os.path.join(parent_dir,"data/examples/5ZOA.gro"), "data/5ZOA.json"), 
          (os.path.join(parent_dir,"data/examples/barstar.gro"), "data/barstar.json"), 
          (os.path.join(parent_dir,"data/examples/DMPC_PI.gro"), "data/DMPC_PI.json"), 
          (os.path.join(parent_dir,"data/examples/noriega_AA_CRD_3CAL.gro"), "data/noriega_AA.json"), 
          (os.path.join(parent_dir,"data/examples/noriega_CG_CRD_3CAL.gro"), "data/noriega_CG.json") ]
)
def test_with_param(filepath_test, JSON_filepath_data):
    JSON_filepath = gd.main(filepath_test)
    with open(JSON_filepath, 'r') as f:
        JSON_content_test = json.load(f)

    with open(JSON_filepath_data, 'r') as f:
        JSON_content_data = json.load(f)    
    
    for index in range(len(JSON_content_data["inventory"])):
        assert JSON_content_test["inventory"][index]["id"] == JSON_content_data["inventory"][index]["id"]
        assert JSON_content_test["inventory"][index]["number_of_atoms"] == JSON_content_data["inventory"][index]["number_of_atoms"]
        assert JSON_content_test["inventory"][index]["number_of_molecules"] == JSON_content_data["inventory"][index]["number_of_molecules"]
        assert JSON_content_test["inventory"][index]["residue_names"] == JSON_content_data["inventory"][index]["residue_names"]
        assert JSON_content_test["inventory"][index]["residue_ids"] == JSON_content_data["inventory"][index]["residue_ids"]
        assert JSON_content_test["inventory"][index]["formula_without_h"] == JSON_content_data["inventory"][index]["formula_without_h"]
        assert JSON_content_test["inventory"][index]["is_solvant"] == JSON_content_data["inventory"][index]["is_solvant"]
        assert JSON_content_test["inventory"][index]["is_ion"] == JSON_content_data["inventory"][index]["is_ion"]
        assert JSON_content_test["inventory"][index]["is_lipid"] == JSON_content_data["inventory"][index]["is_lipid"]
        assert JSON_content_test["inventory"][index]["is_protein"] == JSON_content_data["inventory"][index]["is_protein"]
        assert JSON_content_test["inventory"][index]["protein_sequence"] == JSON_content_data["inventory"][index]["protein_sequence"]
        assert JSON_content_test["inventory"][index]["comment"] == JSON_content_data["inventory"][index]["comment"]
    
    assert JSON_content_test["resolution"] == JSON_content_data["resolution"]
    assert JSON_content_test["file_md5sum"] == JSON_content_data["file_md5sum"]