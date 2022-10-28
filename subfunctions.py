"""###########################################################################
#   This file contains subfunctions for Phase 1 of the TAMU MEEN 357 project
#
#   Created by: MEEN 357 Rover Physics Team
#   Last Modified: 19 September 2022
###########################################################################"""

import math
import numpy as np
from scipy.interpolate import interp1d

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
        output_speed = input_speed / get_gear_ratio(speed_reducer)
        output_torque = input_torque * get_gear_ratio(speed_reducer)

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
    if (type(v) != int) and (type(v) != float) and (not isinstance(v, np.ndarray)):
        raise TypeError('1st argument \'v\' must be a scalar or a numpy array')
    # make v a numpy array if it's a scalar
    if not isinstance(v, np.ndarray):
        v = np.array([v])
    # check that the vector is 1D
    if len(np.shape(v)) != 1:
        raise TypeError('1st argument \'v\' must be a scalar or a vector. Matricies are not allowed.')

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
    if (type(t) != int) and (type(t) != float):
        raise Exception ("argument one must be a scalar")

    # check that y is a np array
    if not isinstance(y, np.ndarray):
        raise Exception("argument two must be a numpy array")
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
    
    terrain_angle = float(alpha_fun(y[1]))
    
    #Veloctiy is first element in array
    velocity = float(y[0])
    
    #calculate accelertation based off a= (1/mass) * Fnet
    a = (1/get_mass(rover)) * F_net(motorW(velocity, rover), terrain_angle, rover, planet, experiment['Crr'])
    
    # 1st element is calculated a, 2nd element is velocity
    dydt[0] = a
    dydt[1] = y[1]
    
    
    return dydt

def mechpower(v, rover):
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
    if (type(v) != int) and (type(v) != float) and (not isinstance(v, np.ndarray)):
        raise TypeError("1st argument must be a scalar or vector")
    
    #check if rover is dict
    if (type(rover) != dict):
        raise TypeError("2nd argument must be a dictionary")
        
    torque = tau_dcmotor(motorW(v, rover), rover['wheel_assembly']['motor'])
    w = motorW(v, rover)
    
    if (type(v) == int) or (type(v) == float):
        
        #find power if given single scalar
        k = 0
        p = torque[k]*w[k]
        
    else:
        
        # check that the vector is 1D
        if len(np.shape(v)) != 1:
            raise TypeError('1st argument \'v\' must be a vector. Matricies are not allowed.')   
        
        # power = torque * omega
        p = np.zeros(len(v), dtype=float)

        #find power for each velocity given
        for k in range(len(v)):
            p[k] = torque[k] * w[k]
            

    return p


def battenergy(t, v, rover): 
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

    # Input Validation

    # check that the time and velocity arguments are numpy arrays
    if not isinstance(t, np.ndarray) or not isinstance(v, np.ndarray):
        raise TypeError('1st and 2nd arguments must be numpy arrays')
    # check that the time and velocity arguments are 1D
    if len(np.shape(t)) != 1 or len(np.shape(v)) != 1:
        raise TypeError('1st and 2nd arguments must be vectors. Matricies are not allowed.')
    # check that the time and velocity arguments are the same size
    if len(t) != len(v):
        raise TypeError('1st and 2nd arguments must be the same size')
    # check that the rover argument is a dict
    if type(rover) != dict:
        raise TypeError('3rd argument must be a dict')
    # check that the rover dictionary contains the required keys
    try:
        motor = rover['wheel_assembly']['motor']
    except KeyError as e:
        raise KeyError(f'Invalid rover dictionary, could not find key: {e}')

    # Main Code

    effcy_tau = rover['wheel_assembly']['motor']['effcy_tau']
    effcy = rover['wheel_assembly']['motor']['effcy']
    effcy_fun = interp1d(effcy_tau, effcy, kind='cubic') # interpolate efficiency data

    # compute the mechanical power output of each motor
    p = mechpower(v, rover)
    print(p)

    # compute the electrical power input to each motor
    omega = motorW(v, rover)
    tau = tau_dcmotor(omega, rover['wheel_assembly']['motor'])
    print(tau)
    
    motor_effcy = effcy_fun(tau)
    e = p / motor_effcy 
    total_e = 6 * e 

    # compute the total electrical energy consumed from the battery pack
    E = np.trapz(total_e, t)

    return E