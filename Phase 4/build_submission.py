"""this script will take the first design from a JSON file and build a submission"""

import json
import os
import sys
from sampling_tools import assemble_design
from subfunctions_Phase4 import redefine_edl_system
from rich import print
import pickle

def build_submission(team_name, team_number, input_path, output_path):
    """main function"""
    with open(input_path, 'r') as input_file:
        designs = json.load(input_file)
    design = designs[0]['design']

    edl = assemble_design(design)
    # just to make sure the system is defined correctly, we will use the same function as in the provided example
    edl = redefine_edl_system(edl)  

    edl['team_name'] = team_name
    edl['team_number'] = team_number

    with open(output_path, 'wb') as handle:
        pickle.dump(edl, handle, protocol=pickle.HIGHEST_PROTOCOL)


if __name__ == '__main__':
    TEAM_NAME = 'Team 1'
    TEAM_NUMBER = 22
    INPUT_PATH = 'Phase 4/optimization_data/multi_terrain_optimization/best_designs.json'
    OUTPUT_PATH = f'Phase 4/submissions/challenge_design_team{TEAM_NUMBER}.pickle'
    build_submission(TEAM_NAME, TEAM_NUMBER, INPUT_PATH, OUTPUT_PATH)