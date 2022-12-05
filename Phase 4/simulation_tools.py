"""
The goal of this file is to provide tools for checking the utility of the sample designs for the rover and edl
"""

import numpy as np
import matplotlib.pyplot as plt
import define_experiment as def_ex
import define_opt_problem as def_opt
import subfunctions_Phase4 as sf
import sampling_tools as st
from rich import print
import scipy

def set_edl_initial_conditions(edl_system):
    edl_system['altitude'] = 11000    # [m] initial altitude
    edl_system['velocity'] = -578     # [m/s] initial velocity
    edl_system['parachute']['deployed'] = True   # our parachute is open
    edl_system['parachute']['ejected'] = False   # and still attached
    edl_system['rover']['on_ground'] = False # the rover has not yet landed
    return edl_system

def check_landing_success(T, Y, edl_system):
    """
    Checks if the landing was successful. Calculates a penalty if it was not.

    Parameters:
        T (array): time array from the simulation. This holdes all the time steps.
        Y (array): state array from the simulation. This holds all the state variables at each time step.
        edl_system (dict): edl system dictionary

    Returns:
        landing_success (bool): True if the landing was successful, False otherwise
        landing_failures (list): list of strings describing the failures
        landing_penalty (float): penalty for the landing
    """
    end_time = T[-1]
    end_velocity = Y[0, -1]
    end_altitude = Y[1, -1]
    end_rover_rel_velocity = Y[5, -1]
    end_rover_rel_position = Y[6, -1]
    end_rover_position = end_altitude + end_rover_rel_position
    end_rover_velocity = end_velocity + end_rover_rel_velocity

    success = True
    failures = []
    penalty_deltas = {}  # will store the delta from the constraint for each constraint that is violated
                         # for example, if the rover traveled 10 m less than the minimum distance, then
                         # penalty_deltas['rover did not make it to destination'] = 10
    # check if end_rover_position is within 1 m of the ground
    if abs(end_rover_position) > 1:
        success = False
        reason_string = 'rover did not land'
        failures.append(reason_string)
        penalty_deltas[reason_string] = abs(end_rover_position - 1)
    if end_altitude <= edl_system["sky_crane"]["danger_altitude"]:
        success = False
        reason_string = 'sky crane too close to ground'
        failures.append(reason_string)
        penalty_deltas[reason_string] = abs(end_altitude - edl_system["sky_crane"]["danger_altitude"])
    if abs(end_rover_velocity) >= abs(edl_system["sky_crane"]["danger_speed"]) and not 'rover did not land' in failures:
        success = False
        reason_string = 'sky crane too fast at landing'
        failures.append(reason_string)
        penalty_deltas[reason_string] = abs(end_rover_velocity - edl_system["sky_crane"]["danger_speed"])
    return success, failures, penalty_deltas

def check_rover_success(edl_system, constraints, tolerance=1e-6):
    """
    Checks if the rover was successful. Calculates a penalty if it was not.

    Parameters:
        telemetry (dict): telemetry dictionary
            'Time':                 time array from the simulation. This holdes all the time steps.
            'completion_time':      time it took to complete the mission
            'velocity':             array of velocities for each time step.
            'position':             array of positions for each time step.
            'distance_traveled':    distance traveled by the rover
            'max_velocity':         maximum velocity of the rover
            'average_velocity':     average velocity of the rover
            'power':                array of power used for each time step.
            'battery_energy':       total battery energy used
            'energy_per_distance':  energy per distance traveled
        constraints (dict): constraints dictionary
    """
    # unpack telemetry
    telemetry = edl_system['rover']['telemetry']

    success = True
    failures = []  # will store strings describing the failures
    penalty_deltas = {}  # will store the delta from the constraint for each constraint that is violated
                         # for example, if the rover traveled 10 m less than the minimum distance, then
                         # penalty_deltas['rover did not make it to destination'] = 10

    # check if the rover made it to the destination
    if abs(telemetry['distance_traveled'] - constraints['distance_traveled'][0]) > tolerance:
        success = False
        reason_string = 'rover did not make it to destination'
        failures.append(reason_string)
        penalty_deltas[reason_string] = abs(telemetry['distance_traveled'] - constraints['distance_traveled'][0])
    # check if the rover ran out of battery
    capacity = edl_system['rover']['power_subsys']['battery']['capacity']
    if telemetry['battery_energy'] > capacity:
        success = False
        reason_string = 'rover ran out of battery'
        failures.append(reason_string)
        penalty_deltas[reason_string] = abs(telemetry['battery_energy'] - capacity)

    return success, failures, penalty_deltas


def evaluate_design(design, experiment, end_event, continuous_options, integer_options, category_options, constraints, tmax=1000, verbose=False, full_output=False):
    # initialize
    edl_system = st.assemble_design(design)
    edl_system = set_edl_initial_conditions(edl_system)
    mission_events = sf.define_mission_events()
    planet = sf.define_planet()
    penalty_weights = def_opt.define_penalty_weights()

    if verbose:
        print(f'[purple]Testing design:[/purple]', design)

    # check if the design is legal
    legal, design_failures, design_penalty_deltas = st.check_design_legality(design, continuous_options, integer_options, category_options, constraints)

    # run the simulations
    if verbose:
        print('[purple]Running edl simulation...[/purple]')
    time_edl_run, edl_end_state, edl_system = sf.simulate_edl(edl_system, planet, mission_events, tmax, verbose)
    time_edl = time_edl_run[-1]
    landing_success, landing_failures, edl_penalty_deltas = check_landing_success(time_edl_run, edl_end_state, edl_system)

    if verbose:
        print('\n[purple]Running rover simulation...[/purple]')
    edl_system['rover'] = sf.simulate_rover(edl_system['rover'], planet, experiment, end_event)
    time_rover = edl_system['rover']['telemetry']['completion_time']
    # print('\nRover telemetry:\n', edl_system["rover"]["telemetry"])
    rover_success, rover_failures, rover_penalty_deltas = check_rover_success(edl_system, constraints)
    

    time = time_edl + time_rover 
    # iterate through the penalty deltas and add the penalty to the total penalty
    design_penalty = sum([design_penalty_deltas[key] * penalty_weights[key] for key in design_penalty_deltas.keys()])
    landing_penalty = sum([edl_penalty_deltas[key] * penalty_weights[key] for key in edl_penalty_deltas.keys()])
    rover_penalty = sum([rover_penalty_deltas[key] * penalty_weights[key] for key in rover_penalty_deltas.keys()])

    total_penalty = design_penalty + landing_penalty + rover_penalty

    if verbose:
        if legal:
            print(f'Design was legal.')
        else:
            print(f'Design was illegal. Failures: {design_failures}')
            print(f'Design penalty deltas: {design_penalty_deltas}')
            print(f'Applied a design penalty of {design_penalty} s.')
        if landing_success:
            print(f'Landing was successful.')
        else:
            print(f'\nLanding was unsuccessful. Failures: {landing_failures}')
            print(f'Landing penalty deltas: {edl_penalty_deltas}')
            print(f'Applied a landing penalty of {landing_penalty} s')
        if rover_success:
            print(f'Rover mission was successful.')
        else:
            print(f'\nRover mission was unsuccessful. Failures: {rover_failures}')
            print(f'Rover mission penalty deltas: {rover_penalty_deltas}')
            print(f'Applied a rover penalty of {rover_penalty} s.')

        if total_penalty == 0:
            print(f'Completed the mission in {time} s.')
        else:
            print(f'Completed the mission in {time} s. Applied a penalty of {total_penalty} s.')
        print(f'\n[dark_blue]Time Result:[/dark_blue] {time}')

    if not full_output:
        return time + total_penalty
    else:
        response = {
            'time': time,
            'success': legal and landing_success and rover_success,
            'legal': legal,
            'landing_success': landing_success,
            'rover_success': rover_success,
            'design_failures': design_failures,
            'landing_failures': landing_failures,
            'rover_failures': rover_failures,
            'design_penalty_deltas': design_penalty_deltas,
            'landing_penalty_deltas': edl_penalty_deltas,
            'rover_penalty_deltas': rover_penalty_deltas,
            'design_penalty': design_penalty,
            'landing_penalty': landing_penalty,
            'rover_penalty': rover_penalty,
            'total_penalty': total_penalty,
            'time_edl': time_edl,
            'time_rover': time_rover,
        }
        return response

def simplified_evaluator_factory(experiment, end_event, continuous_options, integer_options, category_options, constraints, tmax=1000, verbose=False):
    """
    Creates a function to evaluate a design based on just the design num

    Parameters:
        experiment (dict): experiment dictionary
        end_event (dict): end event dictionary
        constraints (dict): constraints dictionary
        tmax (int): maximum time to run the simulation
        verbose (bool): whether to print verbose output

    Returns:
        time (float): time to complete the mission
                        (accounting for penalties)
    """
    def evaluator(design_nums):
        design = st.nums_to_design(design_nums, continuous_options, integer_options, category_options)
        return evaluate_design(design, experiment, end_event, continuous_options, integer_options, category_options, constraints, tmax, verbose)
    return evaluator

def multi_terrain_evaluator_factory(experiment_holder: st.ExperimentHolder, continuous_options, integer_options, category_options, constraints, tmax=1000, verbose=False, full_output=False):
    def evaluator(design_nums):
        design = st.nums_to_design(design_nums, continuous_options, integer_options, category_options)
        total_time = 0
        for experiment, end_event in experiment_holder.experiments:
            time = evaluate_design(design, experiment, end_event, continuous_options, integer_options, category_options, constraints, tmax, verbose, full_output)
            total_time += time
        return total_time
    return evaluator


if __name__ == '__main__':
    # define the experiment
    experiment, end_event = def_ex.experiment1()

    # define the optimization problem
    constraints = def_opt.define_constraints()
    continuous_options = def_opt.define_continuous_options()
    integer_options = def_opt.define_integer_options()
    category_options = def_opt.define_category_options()

    # define the edl
    # designer = st.design_generator(continuous_options, integer_options, category_options, constraints)
    # nums, design = next(designer)

    design = {
        'parachute_diameter': 14.009163985800384,
        'fuel_mass': 181.31447016782022,
        'wheel_radius': 0.688371143604525,
        'gear_diameter': 0.05044328094806457,
        'chassis_mass': 264.63915675195926,
        'num_bat_modules': 15,
        'motor_type': 'speed',
        'battery_type': 'NiMH',
        'chassis_material': 'steel'
    }

    design_nums = [
            0.0003903967023776711,
            0.9178775026465614,
            0.9984545636291741,
            0.00022237290218479755,
            0.013363873810114113,
            0.45560638412083454,
            0.7536689693441452,
            0.18076712121647265,
            0.11787956316768466
        ]

    # evaluate the edl
    # time = evaluate_design(design, experiment, end_event,continuous_options, integer_options, category_options, constraints, verbose=True)
    experiment_holder = st.ExperimentHolder()
    evaluator = multi_terrain_evaluator_factory(experiment_holder, continuous_options, integer_options, category_options, constraints, tmax = 1000, verbose=True)
    time = evaluator(design_nums)
    print(time)