"""
This file will store information about the optimization problem
Bounds, constraints, etc.
"""

import numpy as np

def define_category_options():
    CATEGORY_OPTIONS = {
        'motor_type': ['base', 'base_he', 'torque', 'torque_he', 'speed', 'speed_he'],
        'battery_type': ['LiFePO4', 'NiMH', 'NiCD', 'PbAcid-1', 'PbAcid-2', 'PbAcid-3'],
        'chassis_material': ['steel', 'magnesium', 'carbon'],
    }
    return CATEGORY_OPTIONS

def define_constraints():  
    CONSTRAINTS = {
        'rover_landing_velocity': (-np.inf, -1),
        'distance_traveled': (1000, np.inf),
        'strength': (40000, np.inf),
        'cost': (0, 7.2e6),
    }
    return CONSTRAINTS

def define_continuous_options():
    CONTINUOUS_OPTIONS = {
        'parachute_diameter': (14, 19),  # m
        'fuel_mass': (100, 290),  # kg
        'wheel_radius': (0.2, 0.7),  # m
        'gear_diameter': (0.05, 0.12),  # m
        'chassis_mass': (250, 800),  # kg
    }
    return CONTINUOUS_OPTIONS

def define_integer_options():
    INTEGER_OPTIONS = {
        'num_bat_modules': (1, 40),  # number of battery modules  !!! this one is self imposed by us the students
                                    # I need this to make num_to_design work
    }
    return INTEGER_OPTIONS

def define_penalty_weights():
    PENALTY_WEIGHTS = {
        'rover did not make it to destination': 10,
        'rover ran out of battery': 1e-2,
        'rover did not land': 100,
        'sky crane too close to ground': 100,
        'sky crane too fast at landing': 100,
        'strength is too low': 100,
        'cost is too high': 100,
    }
    return PENALTY_WEIGHTS