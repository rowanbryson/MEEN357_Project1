"""
The goal of this file is to provide tools for creating sample designs for the rover and edl
to be used in the optimization process.
"""

import numpy as np
import matplotlib.pyplot as plt
import define_experiment as def_ex
import subfunctions_Phase4 as sf
from rich import print
import scipy
from scipy.stats.qmc import LatinHypercube

DISCRETE_OPTIONS = {
    'motor': ['base', 'base_he', 'torque', 'torque_he', 'speed', 'speed_he'],
    'battery': ['LiFePO4', 'NiMH', 'NiCD', 'PbAcid-1', 'PbAcid-2', 'PbAcid-3'],
    'chassis': ['steel', 'magnesium', 'carbon'],
}

CONSTRAINTS = {
    'rover_landing_velocity': (-np.inf, -1),
    'distance_traveled': (1000, np.inf),
    'strength': (40000, np.inf),
    'cost': (0, 7.2e6),
}

BOUNDS = {
    'parachute_diameter': (14, 19),  # m
    'fuel_mass': (100, 290),  # kg
    'wheel_radius': (0.2, 0.7),  # m
    'gear_diameter': (0.05, 0.12),  # m
    'chassis_mass': (250, 800),  # kg
    'num_bat_modules': (1, 40),  # number of battery modules  !!! this one is self imposed by us the students
                                 # I need this to make num_to_design work
}

def assemble_design(design: dict) -> dict:
    """
    Assemble the edl from the design dictionary.

    Parameters:
        design (dict): dictionary containing the design parameters
            parachute_diameter (float): diameter of the parachute [m]
            fuel_mass (float): mass of the fuel [kg]
            wheel_radius (float): radius of the wheel [m]
            gear_diameter (float): diameter of the gear [m]
            chassis_mass (float): mass of the chassis [kg]
            motor_type (str): type of motor
            battery_type (str): type of battery
            num_bat_modules (int): number of battery modules
            chassis_material (str): material of the chassis
    """
    edl = sf.define_edl_system()
    edl = sf.define_chassis(edl, design['chassis_material'])
    edl = sf.define_batt_pack(edl, design['battery_type'], design['num_bat_modules'])
    edl = sf.define_motor(edl, design['motor_type'])
    edl['parachute']['diameter'] = design['parachute_diameter']
    edl['rocket']['fuel_mass'] = design['fuel_mass']
    edl['rover']['wheel_assembly']['wheel']['radius'] = design['wheel_radius']
    edl['rover']['wheel_assembly']['speed_reducer']['diam_gear'] = design['gear_diameter']
    edl['rover']['chassis']['mass'] = design['chassis_mass']
    return edl


def check_design_legality(design, bounds):
    """
    Check if the design is legal.
    """
    reason = 'all good'
    for key, value in design.items():
        if key in bounds:
            if value < bounds[key][0] or value > bounds[key][1]:
                reason = f'{key} out of bounds'
                return False
    return True, reason

def nums_to_design(nums, discrete_options, bounds):
    """
    makes a design dictionary based on a tuple of floats between 0 and 1
    """
    design = {}
    for i, key in enumerate(bounds.keys()):
        design[key] = bounds[key][0] + nums[i] * (bounds[key][1] - bounds[key][0])
    for i, key in enumerate(discrete_options.keys()):
        num_choices = len(discrete_options[key])
        design[key] = discrete_options[key][int(nums[len(bounds) + i] * num_choices)]
    return design

def generate_designs(num_samples, discrete_options, bounds):
    # use latin hypercube sampling to create a sample of designs
    # the sample is a list of tuples of floats between 0 and 1
    # the first len(bounds) floats are for the continuous variables
    # the last len(discrete_options) floats are for the discrete variables
    sample = LatinHypercube(d=len(bounds) + len(discrete_options), seed=12345).random(num_samples)
    designs = []
    for nums in sample:
        design = nums_to_design(nums, discrete_options, bounds)
        designs.append(design)
    return designs



if __name__ == '__main__':
    designs = generate_designs(100, DISCRETE_OPTIONS, BOUNDS)
    all_legal = True
    for design in designs:
        legal, reason = check_design_legality(design, BOUNDS)
        if not legal:
            all_legal = False
            print(design)
            print(reason)
    if all_legal:
        print('all designs are legal')