"""
The goal of this file is to provide tools for creating sample designs for the rover and edl
to be used in the optimization process.
"""

import numpy as np
import matplotlib.pyplot as plt
import define_experiment as def_ex
import define_opt_problem as def_opt
import subfunctions_Phase4 as sf
from rich import print
import scipy
from scipy.stats.qmc import LatinHypercube

class ExperimentHolder():

    __atr__ = ['alpha_dist', 'distribution_type', 'alpha_bounds', 'terrain_samples', 'experiments']

    def __init__(self, alpha_dist = None, distribution_type = 'uniform', alpha_bounds = (-5, 20)):
        self.alpha_dist = alpha_dist or np.array([0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000])
        self.distribution_type = distribution_type
        self.alpha_bounds = alpha_bounds

        self.update_experiments(5)

    def update_experiments(self, num_samples):
        """
        Updates the terrain samples to be used in the optimization process.

        Parameters:
            num_samples (int): number of terrain samples to generate
        """
        terrain_samples = [self.generate_terrain() for _ in range(num_samples)]

        experiments = []
        for terrain in terrain_samples:
            base_experiment, base_end_event = def_ex.experiment1()
            base_experiment['alpha_deg'] = terrain['alpha_deg']
            base_experiment['alpha_dist'] = terrain['alpha_dist']
            experiments.append((base_experiment, base_end_event))
        
        self.experiments = experiments


    def generate_terrain(self):
        """
        Generates a terrain map for the rover to traverse.

        Parameters:
            alpha_dist (list): list of floats representing the distances between the alpha values
            distribution_type (str): type of distribution to use for the alpha values
                uniform: uniform distribution
                normal: normal distribution
            alpha_bounds (tuple): tuple of floats representing the bounds for the alpha values

        Returns:
            terrain (dict): dictionary containing the terrain parameters
        """
        alpha_bounds = self.alpha_bounds
        distribution_type = self.distribution_type
        alpha_dist = self.alpha_dist

        if alpha_dist is None:
            alpha_dist = np.linspace(0, 1000, 11)
        if distribution_type == 'uniform':
            alpha_deg = np.random.uniform(alpha_bounds[0], alpha_bounds[1], len(alpha_dist))
            return {'alpha_deg': alpha_deg, 'alpha_dist': alpha_dist}
        elif distribution_type == 'normal':
            alpha_deg = np.random.normal(np.mean(alpha_bounds), np.std(alpha_bounds), len(alpha_dist))
            return {'alpha_deg': alpha_deg, 'alpha_dist': alpha_dist}
        else:
            raise ValueError(f'invalid distribution_type: {distribution_type}')

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
    
    Returns:
        edl (dict): dictionary containing the edl parameters
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


def check_design_legality(design, continuous_options, integer_options, category_options, constraints):
    """
    Check if the design is legal.

    Parameters:
        design (dict): dictionary containing the design parameters
        continuous_options (dict): dictionary of bounds for continuous variables
        discrete_options (dict): dictionary of discrete options
        category_options (dict): dictionary of category options
        constraints (dict): dictionary of constraints

    Returns:
        legal (bool): True if the design is legal, False otherwise
        reason (str): reason why the design is illegal, or None if the design is legal
        penalty_deltas (dict): dictionary of penalty deltas for each constraint
    """
    failures = []
    for key in continuous_options:
        if design[key] < continuous_options[key][0] or design[key] > continuous_options[key][1]:
            failures.append('continuous variable out of bounds')
    for key in integer_options:
        if design[key] < integer_options[key][0] or design[key] > integer_options[key][1]:
            failures.append('integer variable out of bounds')
    for key in category_options:
        if design[key] not in category_options[key]:
            failures.append('category variable out of bounds')

    try:
        edl = assemble_design(design)
    except Exception as e:
        failures.append(f'assemble_design failed: {e}')

    penalty_deltas = {}
    strength = sf.get_cost_edl(edl)
    if strength < constraints['strength'][0] or strength > constraints['strength'][1]:
        reason_string = 'strength is too low'
        failures.append(reason_string)
        penalty_deltas[reason_string] = abs(constraints['strength'][0] - strength)
    cost = sf.get_cost_edl(edl)
    if cost < constraints['cost'][0] or cost > constraints['cost'][1]:
        reason_string = 'cost is too high'
        failures.append(reason_string)
        penalty_deltas[reason_string] = abs(constraints['cost'][1] - cost)

    legal = len(failures) == 0
    return legal, failures, penalty_deltas

def nums_to_design(nums, continuous_options, integer_options, category_options):
    """
    makes a design dictionary based on a tuple of floats between 0 and 1

    we have 9 degrees of freedom total, with 6 continuous and 3 discrete
    the first 6 floats are for the continuous variables
    the last 3 floats are for the discrete variables

    Parameters:
        nums (tuple): tuple of floats between 0 and 1
            len(nums) must equal len(bounds) + len(discrete_options) + len(category_options)
        discrete_options (dict): dictionary of discrete options
        bounds (dict): dictionary of bounds for continuous variables

    Returns:
        design (dict): dictionary containing the design parameters
    """
    design = {}
    for i, key in enumerate(continuous_options):
        design[key] = continuous_options[key][0] + nums[i] * (continuous_options[key][1] - continuous_options[key][0])
    for i, key in enumerate(integer_options):
        design[key] = int(integer_options[key][0] + nums[i + len(continuous_options)] * (integer_options[key][1] - integer_options[key][0]))
    for i, key in enumerate(category_options):
        design[key] = category_options[key][int(nums[i + len(continuous_options) + len(integer_options)] * len(category_options[key]))]
    return design

def design_generator(continuous_options, integer_options, category_options, constraints):
    """
    Infinite generator of legal designs.

    Parameters:
        continuous_options (dict): dictionary of bounds for continuous variables
        discrete_options (dict): dictionary of discrete options
        category_options (dict): dictionary of category options
        constraints (dict): dictionary of constraints

    Yields:
        design_nums (tuple): tuple of floats between 0 and 1 representing the design
        design (dict): dictionary containing the design parameters
    """
    # use latin hypercube sampling to create a sample of designs
    # the sample is a list of tuples of floats between 0 and 1

    degrees_of_freedom = len(continuous_options) + len(integer_options) + len(category_options)
    hypercube = LatinHypercube(d=degrees_of_freedom)
    while True:
        design_nums = hypercube.random()[0]
        design = nums_to_design(design_nums, continuous_options, integer_options, category_options)
        legal, _, _ = check_design_legality(design, continuous_options, integer_options, category_options, constraints)
        if legal:
            yield design_nums, design

def terrain_generator(alpha_dist = None, distribution_type = 'uniform', alpha_bounds = (-5, 20)):
    """
    Infinite generator of terrain.

    Parameters:
        alpha_dist (list): list of floats representing the distances between the alpha values
        distribution_type (str): type of distribution to use for the alpha values
            uniform: uniform distribution
            normal: normal distribution
        alpha_bounds (tuple): tuple of floats representing the bounds for the alpha values

    Yields:
        terrain (dict): dictionary containing the terrain parameters
    """
    if alpha_dist is None:
        alpha_dist = np.linspace(0, 1000, 11)

    if distribution_type == 'uniform':
        while True:
            alpha_deg = np.random.uniform(alpha_bounds[0], alpha_bounds[1], len(alpha_dist))
            yield {'alpha_deg': alpha_deg, 'alpha_dist': alpha_dist}
    elif distribution_type == 'normal':
        while True:
            alpha_deg = np.random.normal(np.mean(alpha_bounds), np.std(alpha_bounds), len(alpha_dist))
            yield {'alpha_deg': alpha_deg, 'alpha_dist': alpha_dist}
    else:
        raise ValueError(f'invalid distribution_type: {distribution_type}')


if __name__ == '__main__':
    # initialize the experiment
    CONTINUOUS_OPTIONS = def_opt.define_continuous_options()
    INTEGER_OPTIONS = def_opt.define_integer_options()
    CATEGORY_OPTIONS = def_opt.define_category_options()
    CONSTRAINTS = def_opt.define_constraints()

    designer = design_generator(CONTINUOUS_OPTIONS, INTEGER_OPTIONS, CATEGORY_OPTIONS, CONSTRAINTS)
    for i in range(10):
        design = next(designer)
        print(design)