"""###########################################################################
#   This file contains subfunctions for Phase 1 of the TAMU MEEN 357 project
#
#   Created by: MEEN 357 Rover Physics Team
#   Last Modified: 19 September 2022
###########################################################################"""

import math
from tkinter import N
from tkinter.ttk import Style
import numpy as np
from scipy.interpolate import interp1d
from functools import partial
from scipy.integrate import solve_ivp, simpson, cumulative_trapezoid
from matplotlib import pyplot as plt
from typing import Union

def get_mass(rover):
    """
    Inputs:  rover:  dict      Data structure containing rover parameters
    
    Outputs:     m:  scalar    Rover mass [kg].
    """
    
    # Check that the input is a dict
    if type(rover) != dict:
        raise Exception('Input must be a dict')
    
    # add up mass of chassis, power subsystem, science payload, 
    # and components from all six wheel assemblies
    m = rover['chassis']['mass'] \
        + rover['power_subsys']['mass'] \
        + rover['science_payload']['mass'] \
        + 6*rover['wheel_assembly']['motor']['mass'] \
        + 6*rover['wheel_assembly']['speed_reducer']['mass'] \
        + 6*rover['wheel_assembly']['wheel']['mass'] \
    
    return m


def get_gear_ratio(speed_reducer):
    """
    !!! examples:
        wheel_W = motor_W * get_gear_ratio(speed_reducer)
        motor_W = wheel_W / get_gear_ratio(speed_reducer)
        wheel_torque = motor_torque / get_gear_ratio(speed_reducer)
        motor_torque = wheel_torque * get_gear_ratio(speed_reducer)

    Inputs:  speed_reducer:  dict      Data dictionary specifying speed
                                        reducer parameters
    Outputs:            Ng:  scalar    Speed ratio from input pinion shaft
                                        to output gear shaft. Unitless.
    """
    
    
    # Check that the input is a dict
    if type(speed_reducer) != dict:
        raise Exception('Input must be a dict')
    
    # Check 'type' field (not case sensitive)
    if speed_reducer['type'].lower() != 'reverted':
        raise Exception('The speed reducer type is not recognized.')
    
    # Main code
    d1 = speed_reducer['diam_pinion']
    d2 = speed_reducer['diam_gear']
    
    Ng = (d2/d1)**2
    
    return Ng


def tau_dcmotor(omega, motor):
    """
    Inputs:  omega:  numpy array      Motor shaft speed [rad/s]
             motor:  dict             Data dictionary specifying motor parameters
    Outputs:   tau:  numpy array      Torque at motor shaft [Nm].  Return argument
                                      is same size as first input argument.
    """

    
    # Check that the first input is a scalar or a vector
    if (type(omega) != int) and (type(omega) != float) and (not isinstance(omega, np.ndarray)):
        raise Exception('First input must be a scalar or a vector. If input is a vector, it should be defined as a numpy array.')
    elif not isinstance(omega, np.ndarray):
        omega = np.array([omega],dtype=float) # make the scalar a numpy array
    elif len(np.shape(omega)) != 1:
        raise Exception('First input must be a scalar or a vector. Matrices are not allowed.')

    # Check that the second input is a dict
    if type(motor) != dict:
        raise Exception('Second input must be a dict')
        
    # Main code
    tau_s    = motor['torque_stall']
    tau_nl   = motor['torque_noload']
    omega_nl = motor['speed_noload']
    
    # initialize
    tau = np.zeros(len(omega),dtype = float)
    for ii in range(len(omega)):
        if omega[ii] >= 0 and omega[ii] <= omega_nl:
            tau[ii] = tau_s - (tau_s-tau_nl)/omega_nl *omega[ii]
        elif omega[ii] < 0:
            tau[ii] = tau_s
        elif omega[ii] > omega_nl:
            tau[ii] = 0
        
    return tau
    
    


def F_rolling(omega, terrain_angle, rover, planet, Crr):
    """
    Inputs:           omega:  numpy array     Motor shaft speed [rad/s]
              terrain_angle:  numpy array     Array of terrain angles [deg]
                      rover:  dict            Data structure specifying rover 
                                              parameters
                    planet:  dict            Data dictionary specifying planetary 
                                              parameters
                        Crr:  scalar          Value of rolling resistance coefficient
                                              [-]
    
    Outputs:           Frr:  numpy array     Array of forces [N]
    """
    
    # Check that the first input is a scalar or a vector
    if (type(omega) != int) and (type(omega) != float) and (not isinstance(omega, np.ndarray)):
        raise Exception('First input must be a scalar or a vector. If input is a vector, it should be defined as a numpy array.')
    elif not isinstance(omega, np.ndarray):
        omega = np.array([omega],dtype=float) # make the scalar a numpy array
    elif len(np.shape(omega)) != 1:
        raise Exception('First input must be a scalar or a vector. Matrices are not allowed.')
        
    # Check that the second input is a scalar or a vector
    if (type(terrain_angle) != int) and (type(terrain_angle) != float) and (not isinstance(terrain_angle, np.ndarray)):
        raise Exception('Second input must be a scalar or a vector. If input is a vector, it should be defined as a numpy array.')
    elif not isinstance(terrain_angle, np.ndarray):
        terrain_angle = np.array([terrain_angle],dtype=float) # make the scalar a numpy array
    elif len(np.shape(terrain_angle)) != 1:
        raise Exception('Second input must be a scalar or a vector. Matrices are not allowed.')
        
    # Check that the first two inputs are of the same size
    if len(omega) != len(terrain_angle):
        raise Exception('First two inputs must be the same size')
    
    # Check that values of the second input are within the feasible range  
    if max([abs(x) for x in terrain_angle]) > 75:    
        raise Exception('All elements of the second input must be between -75 degrees and +75 degrees')
        
    # Check that the third input is a dict
    if type(rover) != dict:
        raise Exception('Third input must be a dict')
        
    # Check that the fourth input is a dict
    if type(planet) != dict:
        raise Exception('Fourth input must be a dict')
        
    # Check that the fifth input is a scalar and positive
    if (type(Crr) != int) and (type(Crr) != float):
        raise Exception('Fifth input must be a scalar')
    if Crr <= 0:
        raise Exception('Fifth input must be a positive number')
        
    # Main Code
    m = get_mass(rover)
    g = planet['g']
    r = rover['wheel_assembly']['wheel']['radius']
    Ng = get_gear_ratio(rover['wheel_assembly']['speed_reducer'])
    
    v_rover = r*omega/Ng
    
    Fn = np.array([m*g*math.cos(math.radians(x)) for x in terrain_angle],dtype=float) # normal force
    Frr_simple = -Crr*Fn # simple rolling resistance
    
    Frr = np.array([math.erf(40*v_rover[ii]) * Frr_simple[ii] for ii in range(len(v_rover))], dtype = float)
    
    return Frr


def F_gravity(terrain_angle, rover, planet):
    """
    Inputs:  terrain_angle:  numpy array   Array of terrain angles [deg]
                     rover:  dict          Data structure specifying rover 
                                            parameters
                    planet:  dict          Data dictionary specifying planetary 
                                            parameters
    
    Outputs:           Fgt:  numpy array   Array of forces [N]
    """
    
    # Check that the first input is a scalar or a vector
    if (type(terrain_angle) != int) and (type(terrain_angle) != float) and (not isinstance(terrain_angle, np.ndarray)):
        raise Exception('First input must be a scalar or a vector. If input is a vector, it should be defined as a numpy array.')
    elif not isinstance(terrain_angle, np.ndarray):
        terrain_angle = np.array([terrain_angle],dtype=float) # make the scalar a numpy array
    elif len(np.shape(terrain_angle)) != 1:
        raise Exception('First input must be a scalar or a vector. Matrices are not allowed.')
        
    # Check that values of the first input are within the feasible range  
    if max([abs(x) for x in terrain_angle]) > 75:    
        raise Exception('All elements of the first input must be between -75 degrees and +75 degrees')

    # Check that the second input is a dict
    if type(rover) != dict:
        raise Exception('Second input must be a dict')
    
    # Check that the third input is a dict
    if type(planet) != dict:
        raise Exception('Third input must be a dict')
        
    # Main Code
    m = get_mass(rover)
    g = planet['g']
    
    Fgt = np.array([-m*g*math.sin(math.radians(x)) for x in terrain_angle], dtype = float)
        
    return Fgt


def F_drive(omega, rover):
    """
    Inputs:  omega:  numpy array   Array of motor shaft speeds [rad/s]
             rover:  dict          Data dictionary specifying rover parameters
    
    Outputs:    Fd:  numpy array   Array of drive forces [N]
    """
    
    
    # Check that the first input is a scalar or a vector
    if (type(omega) != int) and (type(omega) != float) and (not isinstance(omega, np.ndarray)):
        raise Exception('First input must be a scalar or a vector. If input is a vector, it should be defined as a numpy array.')
    elif not isinstance(omega, np.ndarray):
        omega = np.array([omega],dtype=float) # make the scalar a numpy array
    elif len(np.shape(omega)) != 1:
        raise Exception('First input must be a scalar or a vector. Matrices are not allowed.')

    # Check that the second input is a dict
    if type(rover) != dict:
        raise Exception('Second input must be a dict')
    
    # Main code
    Ng = get_gear_ratio(rover['wheel_assembly']['speed_reducer'])
    
    tau = tau_dcmotor(omega, rover['wheel_assembly']['motor'])
    tau_out = tau*Ng
    
    r = rover['wheel_assembly']['wheel']['radius']
    
    # Drive force for one wheel
    Fd_wheel = tau_out/r 
    
    # Drive force for all six wheels
    Fd = 6*Fd_wheel
    
    return Fd


def F_net(omega, terrain_angle, rover, planet, Crr):
    """
    Inputs:           omega:  numpy array     Motor shaft speed [rad/s]
              terrain_angle:  numpy array     Array of terrain angles [deg]
                      rover:  dict     Data structure specifying rover 
                                      parameters
                     planet:  dict     Data dictionary specifying planetary 
                                      parameters
                        Crr:  scalar   Value of rolling resistance coefficient
                                      [-]
    
    Outputs:           Fnet:  numpy array     Array of forces [N]
    """
    
    # Check that the first input is a scalar or a vector
    if (type(omega) != int) and (type(omega) != float) and (not isinstance(omega, np.ndarray)):
        raise Exception('First input must be a scalar or a vector. If input is a vector, it should be defined as a numpy array.')
    elif not isinstance(omega, np.ndarray):
        omega = np.array([omega],dtype=float) # make the scalar a numpy array
    elif len(np.shape(omega)) != 1:
        raise Exception('First input must be a scalar or a vector. Matrices are not allowed.')
        
    # Check that the second input is a scalar or a vector
    if (type(terrain_angle) != int) and (type(terrain_angle) != float) and (not isinstance(terrain_angle, np.ndarray)):
        raise Exception('Second input must be a scalar or a vector. If input is a vector, it should be defined as a numpy array.')
    elif not isinstance(terrain_angle, np.ndarray):
        terrain_angle = np.array([terrain_angle],dtype=float) # make the scalar a numpy array
    elif len(np.shape(terrain_angle)) != 1:
        raise Exception('Second input must be a scalar or a vector. Matrices are not allowed.')
        
    # Check that the first two inputs are of the same size
    if len(omega) != len(terrain_angle):
        raise Exception('First two inputs must be the same size')
    
    # Check that values of the second input are within the feasible range  
    if max([abs(x) for x in terrain_angle]) > 75:    
        raise Exception('All elements of the second input must be between -75 degrees and +75 degrees')
        
    # Check that the third input is a dict
    if type(rover) != dict:
        raise Exception('Third input must be a dict')
        
    # Check that the fourth input is a dict
    if type(planet) != dict:
        raise Exception('Fourth input must be a dict')
        
    # Check that the fifth input is a scalar and positive
    if (type(Crr) != int) and (type(Crr) != float):
        raise Exception('Fifth input must be a scalar')
    if Crr <= 0:
        raise Exception('Fifth input must be a positive number')
    
    # Main Code
    Fd = F_drive(omega, rover)
    Frr = F_rolling(omega, terrain_angle, rover, planet, Crr)
    Fg = F_gravity(terrain_angle, rover, planet)
    
    Fnet = Fd + Frr + Fg # signs are handled in individual functions
    
    return Fnet


def end_of_mission_event(end_event):
    """
    Defines an event that terminates the mission simulation. Mission is over
    when rover reaches a certain distance, has moved for a maximum simulation 
    time or has reached a minimum velocity.            
    """
    
    mission_distance = end_event['max_distance']
    mission_max_time = end_event['max_time']
    mission_min_velocity = end_event['min_velocity']
    
    # Assume that y[1] is the distance traveled
    distance_left = lambda t,y: mission_distance - y[1]
    distance_left.terminal = True
    
    time_left = lambda t,y: mission_max_time - t
    time_left.terminal = True
    
    velocity_threshold = lambda t,y: y[0] - mission_min_velocity
    velocity_threshold.terminal = True
    velocity_threshold.direction = -1
    
    # terminal indicates whether any of the conditions can lead to the
    # termination of the ODE solver. In this case all conditions can terminate
    # the simulation independently.
    
    # direction indicates whether the direction along which the different
    # conditions is reached matter or does not matter. In this case, only
    # the direction in which the velocity treshold is arrived at matters
    # (negative)
    
    events = [distance_left, time_left, velocity_threshold]
    
    return events


def motorW(v, rover):
    '''
    Compute the rotational speed of the motor shaft [rad/s] given the translational velocity of the rover and the rover
    dictionary.
    This function should be “vectorized” such that if given a vector of rover velocities it returns a vector the same size
    containing the corresponding motor speeds.

    Inputs
    ------
    v: scalar or 1D numpy array
        Rover translational velocity [m/s]
    rover: dict
        Data structure containing rover parameters

    Outputs
    -------
    w: scalar or 1D numpy array
        Motor speed [rad/s]
        Return argument should match type/size of input
    '''

    # check that the velocity argument is a scalar or numpy array
    if (type(v) != int) and (type(v) != float) and (np.ScalarType) and (not isinstance(v, np.ndarray)):
        raise TypeError('1st argument \'v\' must be a scalar or a numpy array')
    # make v a numpy array if it's a scalar
    # if not isinstance(v, np.ndarray):
    #     v = np.array([v])
    # check that the vector is 1D
    if isinstance(v, np.ndarray) and len(np.shape(v)) != 1:
        raise TypeError(f'1st argument \'v\' must be a scalar or a vector, not {type(v)} Matricies are not allowed. np.shape(v) = {np.shape(v)}')

    try:
        gear_ratio = get_gear_ratio(rover['wheel_assembly']['speed_reducer'])  # ratio of gearbox input shaft speed to output shaft speed
        wheel_w = v / rover['wheel_assembly']['wheel']['radius']  # rotational velocity of the wheel [rad/s]
    except KeyError as e:
        raise KeyError(f'Invalid rover dictionary, could not find key: {e}')

    motor_w = wheel_w * gear_ratio

    return motor_w

def rover_dynamics(t, y, rover, planet, experiment):
    '''
    Compute the derivative of the [velocity, position] vector for the rover given it's current state

    Inputs
    ------
    t: scalar
        Time sample [s]

    y: 1D numpy array
        Two element array of dependent variables
        First element is rover velocity[m/s], second is rover position [m]

    rover: Dict
        Data structure containing rover definition

    planet: Dict
        Data structure containing planet definition

    experiment: Dict
        Data structure containing experiment defintion

    Outputs
    -------
    dydt: 1D numpy array
        Two element array of first derivatives of state vector
        First element is rover acceleration [m/s^2] and second is rover velocity [m/s]  
    '''
    # check that t is a scalar
    if (type(t) != int) and (type(t) != float) and (not isinstance(t, np.ScalarType)):
        raise Exception (f"argument one must be a scalar, not {type(t)}")

    # check that y is a np array
    if not isinstance(y, np.ndarray):
        raise Exception(f"argument two must be a numpy array, not {type(y)}")
    # check that the vector is 1D
    if len(np.shape(y)) != 1:
        raise TypeError('2nd argument \'y\' must be a vector. Matricies are not allowed.')
    #check that y only has two elements
    if len(y) != 2:
        raise Exception("y must only contain two elements, acceleration and velocity")

    # check that rover, planet, and experiment are dicitonaries 
    if (type(rover) != dict) or (type(planet) != dict) or (type(experiment) != dict):
        raise Exception("Third, fourth and fifth arguments must be dictionaries")
        
    # define alpha_fun
    alpha_fun = interp1d(experiment['alpha_dist'], experiment['alpha_deg'], kind = 'cubic', fill_value = 'extrapolate') 
    
    dydt = np.zeros(len(y), dtype=float)
    
    #terrain angle is interpolated position
    terrain_angle = float(alpha_fun(y[1]))
    
    #Veloctiy is first element in array
    velocity = float(y[0])
    
    #calculate accelertation based off a= (1/mass) * Fnet
    a = (1/get_mass(rover)) * F_net(motorW(velocity, rover), terrain_angle, rover, planet, experiment['Crr'])
    
    # 1st element is calculated a, 2nd element is velocity
    dydt[0] = a
    dydt[1] = y[0]
    
    return dydt

def mechpower(v, rover, quick_plot=False):
    '''
    Computes the instantaneous mechanical power output by a single DC motor at each point in a given velocity profile

    Inputs
    ------
    v: 1D numpy array
        Rover velocity data from a simulation [m/s]
    rover: Dict
        Data structure containing rover definition

    Outputs
    -------
    p: 1D numpy array or scalar
        Instantneous power output of a single motor corresponding to each element in v [W]
        Return argument should be the same size as v
    '''
    #check if v is scalar or vector
    if (type(v) != int) and (type(v) != float) and (not isinstance(v, np.ndarray)) and (not isinstance(v, np.ScalarType)):
        raise TypeError(f"1st argument must be a scalar or vector, not {type(v)}")
    
    #check if rover is dict
    if (type(rover) != dict):
        raise TypeError("2nd argument must be a dictionary")
        
    INPUT_TYPE = type(v)
    if not isinstance(v, np.ndarray):
        v = np.array([v])

    torque = tau_dcmotor(motorW(v, rover), rover['wheel_assembly']['motor'])
    w = motorW(v, rover)
    power = torque * w

    # if (type(v) == int) or (type(v) == float) or isinstance(v, np.ScalarType):
        
    #     #find power if given single scalar
    #     k = 0
    #     p = torque[k]*w[k]
        
    # else:
        
    #     # check that the vector is 1D
    #     if len(np.shape(v)) != 1:
    #         raise TypeError('1st argument \'v\' must be a vector. Matricies are not allowed.')   
        
    #     # power = torque * omega
    #     p = np.zeros(len(v), dtype=float)

    #     #find power for each velocity given
    #     for k in range(len(v)):
    #         p[k] = torque[k] * w[k]
            
    if not INPUT_TYPE == np.ndarray:
        power = power[0]

    if quick_plot:
        fig, ax = plt.subplots(4, 1)
        ax[0].plot(v, torque)
        ax[0].set_xlabel('Velocity [m/s]')
        ax[0].set_ylabel('Torque [Nm]')
        ax[1].plot(w, torque)
        ax[1].set_xlabel('Motor Speed [rad/s]')
        ax[1].set_ylabel('Torque [Nm]')
        ax[2].plot(v, w)
        ax[2].set_xlabel('Velocity [m/s]')
        ax[2].set_ylabel('Motor Speed [rad/s]')
        ax[3].plot(v, power)
        ax[3].set_xlabel('Velocity [m/s]')
        ax[3].set_ylabel('Power [W]')
        fig.tight_layout()
        plt.show()

    return power


def battenergy(t, v, rover, quick_plot=False): 
    """
    This function computes the total electrical energy consumed from the rover battery pack over a simulation profile, 
    defined as time-velocity pairs. This function assumes all 6 motors are driven from the same battery pack (i.e., this 
    function accounts for energy consumed by all motors). 

    Inputs:
    t: 1D numpy array
        Time vector [s]
    v: 1D numpy array
        Velocity vector [m/s]
    rover: dict
        Data structure containing rover parameters

    Outputs:
    E: scalar
        Total electrical energy consumed from the battery pack [J]
    """
    # get an interpolated function for the battery efficiency
    motor = rover['wheel_assembly']['motor']
    motor_effcy_func = interp1d(motor['effcy_tau'], motor['effcy'], kind='cubic', fill_value='extrapolate')

    # get the total mechanical power output by the motors at each velocity
    mechanical_power_vals = 6 * mechpower(v, rover)

    if quick_plot:
        fig, ax = plt.subplots()
        ax.plot(motor['effcy_tau'], motor['effcy'], label='Data')
        interp_x = np.linspace(motor['effcy_tau'][0], motor['effcy_tau'][-1], 1000)
        ax.plot(interp_x, motor_effcy_func(interp_x), label='Interpolation')
        ax.set_title('Motor Efficiency vs. Mechanical Power')
        ax.set_xlabel('Mechanical Power [W]')
        ax.set_ylabel('Motor Efficiency')
        ax.legend()
        plt.show()

    # get the battery efficiency at each mechanical power value
    motor_effcy = motor_effcy_func(tau_dcmotor(motorW(v, rover), motor))
    # get the total electrical power consumed by the motors at each velocity
    electrical_power = mechanical_power_vals / motor_effcy

    # integrate the electrical power to get the total electrical energy consumed
    trapazoid_answer = cumulative_trapezoid(electrical_power, t)[-1] # trapanzoidal answer is closer to the value provided by prof
    simpson_answer = simpson(electrical_power, t)

    if quick_plot:
        print(f'Trapazoid answer: {trapazoid_answer}')
        print(f'Simpson answer: {simpson_answer}')
        fig, ax = plt.subplots()
        ax.scatter(t, electrical_power)
        ax.set_title('Electrical Power vs. Time')
        ax.set_xlabel('Time [s]')
        ax.set_ylabel('Electrical Power [W]')
        plt.show()
    
    return simpson_answer

def simulate_rover(rover: dict, planet: dict, experiment: dict, end_event: dict=None, quick_plot: bool=False) -> dict:
    '''
    This function integrates the trajectory of a rover.

    Inputs
    ------
    rover: dict
        Data structure containing the parameters of the rover
    planet: dict
        Data structure containing the planet definition
    experiment: dict
        Data structure containing parameters of the trajectory to be followed by the rover
    end_event: dict (optional)
        Data structure containing the conditions necessary and sufficient to terminate simulation of rover dynamics
        overrides experiment['end_event'] if both are provided
        defaults to experiment['end_event'] if not provided
    
    Outputs
    -------
    rover: dict
        Data structure containing the parameters of the rover, including updated telemetry information.
        telemetry information is stored in the 'telemetry' key of the rover dictionary, with the following keys:

            time: 1D numpy array
                N-element array containing the time history of the rover [s]
            completion_time: scalar
                Time to complete a mission [s]
            velocity: 1D numpy array
                N-element array containing the velocity of the rover as it follows a trajectory [m/s]
            position: 1D numpy array
                N-element array containing the position of the rover as it follows a trajectory [m]
            distance_traveled: scalar
                Total distance traveled by the rover [m]
            max_velocity: scalar
                Maximum velocity of rover along a givien trajectory [m/s]
            average_velocity: scalar
                Average velocity of rover along a given trajectory [m/s]
            power: 1D numpy array
                N-element array containing the instantaneous power outputted by the motor along a trajectory [W]
            battery_energy: scalar
                Total energy to be extracted from the battery to complete trajectory [J]
            energy_per_distance: scalar
                Total energy spent (from battery) per meter traveled [J/m]
    '''

    if end_event is None:
        end_event = experiment['end_event']

    time_span = experiment['time_range']
    y0 = experiment['initial_conditions']

    rover_dynamics_partial = partial(rover_dynamics, rover=rover, planet=planet, experiment=experiment)
    sol = solve_ivp(rover_dynamics_partial, time_span, y0, method='BDF', events=end_of_mission_event(end_event), dense_output=True, max_step=0.1)
    t = sol.t
    y = sol.y

    telemetry = {
        'time': t,
        'completion_time': t[-1],
        'velocity': y[0],
        'position': y[1],
        'distance_traveled': y[1][-1] - y[1][0],
        'max_velocity': np.max(y[0]),
        'average_velocity': (y[1][-1] - y[1][0])/t[-1],
        'power': mechpower(y[0], rover),
        'battery_energy': battenergy(t, y[0], rover),
        'energy_per_distance': battenergy(t, y[0], rover)/(y[1][-1] - y[1][0])
    }
    rover = rover.copy()
    rover['telemetry'] = telemetry
    return rover




    print(sol.message)
    # figure out which event terminal event was triggered
    if sol.status == 1:
        end_event_index = None
        ending_time = t[-1]
        for ii, t_event in enumerate(sol.t_events):
            if len(t_event) > 0:
                if t_event[0] == ending_time:
                    end_event_index = ii
        end_event_reason = list(end_event.keys())[end_event_index]
        print(f'Ending condition: {end_event_reason}')

    if quick_plot:
        import matplotlib.pyplot as plt
        # make 2 by 1 plot of position and velocity
        fig, ax = plt.subplots(2, 1)
        ax[0].plot(t, y[0], color='orange', label='velocity')
        ax[0].set_ylabel('velocity [m/s]')
        ax[0].set_xlabel('time [s]')
        ax[1].plot(t, y[1], color='blue' , label='position')
        ax[1].set_ylabel('position [m]')
        ax[1].set_xlabel('time [s]')
        plt.show()