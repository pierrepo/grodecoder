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
          (os.path.join(parent_dir,"data/examples/noriega_CG_CRD_3CAL.gro"), "data/noriega_CG.json"), 
          (os.path.join(parent_dir,"data/examples/step7_1.gro"), "data/step7_1.json")]
)
def test_with_param(filepath_test, JSON_filepath_data):
    JSON_filepath = gd.main(filepath_test)
    with open(JSON_filepath, 'r') as f:
        JSON_content_test = json.load(f)

    with open(JSON_filepath_data, 'r') as f:
        JSON_content_data = json.load(f)

    assert JSON_content_test["inventory"] == JSON_content_data["inventory"]
