import numpy as np
# import numba as nb
from subfunctions import *
import functools as ft

# @nb.njit()
# def __get_gear_ratio(d1,d2):
#     return (d2/d1) ** 2

# def F_gravity(terrain_angle, rover, planet):
#     """
#     Inputs:  terrain_angle:  numpy array   Array of terrain angles [deg]
#                      rover:  dict          Data structure specifying rover 
#                                             parameters
#                     planet:  dict          Data dictionary specifying planetary 
#                                             parameters
    
#     Outputs:           Fgt:  numpy array   Array of forces [N]
#     """
    
#     # Check that the first input is a scalar or a vector
#     if (type(terrain_angle) != int) and (type(terrain_angle) != float) and (not isinstance(terrain_angle, np.ndarray)):
#         raise Exception('First input must be a scalar or a vector. If input is a vector, it should be defined as a numpy array.')
#     elif not isinstance(terrain_angle, np.ndarray):
#         terrain_angle = np.array([terrain_angle],dtype=float) # make the scalar a numpy array
#     elif len(np.shape(terrain_angle)) != 1:
#         raise Exception('First input must be a scalar or a vector. Matrices are not allowed.')
        
#     # Check that values of the first input are within the feasible range  
#     if max([abs(x) for x in terrain_angle]) > 75:    
#         raise Exception('All elements of the first input must be between -75 degrees and +75 degrees')

#     # Check that the second input is a dict
#     if type(rover) != dict:
#         raise Exception('Second input must be a dict')
    
#     # Check that the third input is a dict
#     if type(planet) != dict:
#         raise Exception('Third input must be a dict')
        
#     # Main Code
#     m = get_mass(rover)
#     g = planet['g']
    
#     Fgt = -m * g * np.sin(terrain_angle * np.pi / 180)

#     # Fgt = np.array([-m*g*math.sin(math.radians(x)) for x in terrain_angle], dtype = float)
        
#     return Fgt

# def get_gear_ratio(speed_reducer):
#     """
#     !!! examples:
#         wheel_W = motor_W * get_gear_ratio(speed_reducer)
#         motor_W = wheel_W / get_gear_ratio(speed_reducer)
#         wheel_torque = motor_torque / get_gear_ratio(speed_reducer)
#         motor_torque = wheel_torque * get_gear_ratio(speed_reducer)

#     Inputs:  speed_reducer:  dict      Data dictionary specifying speed
#                                         reducer parameters
#     Outputs:            Ng:  scalar    Speed ratio from input pinion shaft
#                                         to output gear shaft. Unitless.
#     """
    
    
#     # Check that the input is a dict
#     # if type(speed_reducer) != dict:
#     #     raise Exception('Input must be a dict')
    
#     # # Check 'type' field (not case sensitive)
#     # if speed_reducer['type'].lower() != 'reverted':
#     #     raise Exception('The speed reducer type is not recognized.')
    
#     # Main code
#     d1 = speed_reducer['diam_pinion']
#     d2 = speed_reducer['diam_gear']

#     Ng = __get_gear_ratio(d1,d2)
    
#     return Ng

# @nb.vectorize()
# def __tau_dc_motor(omega, tau_s, tau_nl, omega_nl):
#         if omega >= 0 and omega <= omega_nl:
#             return tau_s - (tau_s-tau_nl)/omega_nl *omega
#         elif omega < 0:
#             return tau_s
#         else:
#             return 0
    

# def tau_dcmotor(omega, motor):
#     """
#     Inputs:  omega:  numpy array      Motor shaft speed [rad/s]
#              motor:  dict             Data dictionary specifying motor parameters
#     Outputs:   tau:  numpy array      Torque at motor shaft [Nm].  Return argument
#                                       is same size as first input argument.
#     """

    
#     # Check that the first input is a scalar or a vector
#     if (type(omega) != int) and (type(omega) != float) and (not isinstance(omega, np.ndarray)):
#         raise Exception('First input must be a scalar or a vector. If input is a vector, it should be defined as a numpy array.')
#     elif not isinstance(omega, np.ndarray):
#         omega = np.array([omega],dtype=float) # make the scalar a numpy array
#     elif len(np.shape(omega)) != 1:
#         raise Exception('First input must be a scalar or a vector. Matrices are not allowed.')

#     # Check that the second input is a dict
#     if type(motor) != dict:
#         raise Exception('Second input must be a dict')
        
    
#     # Main code
#     tau_s    = motor['torque_stall']
#     tau_nl   = motor['torque_noload']
#     omega_nl = motor['speed_noload']
    
#     # call the numba function
#     tau = __tau_dc_motor(omega, tau_s, tau_nl, omega_nl)
        
#     return tau

# @nb.vectorize()
# def __F_rolling(omega, terrain_angle, Crr, m, g, r, Ng):
#     v_rover = r*omega/Ng
    
#     # Fn = np.array([m*g*math.cos(math.radians(x)) for x in terrain_angle],dtype=float) # normal force
#     Fn = m*g*math.cos(math.radians(terrain_angle)) # normal force
#     Frr_simple = -Crr*Fn # simple rolling resistance
    
#     # Frr = np.array([math.erf(40*v_rover[ii]) * Frr_simple[ii] for ii in range(len(v_rover))], dtype = float)
#     Frr = math.erf(40*v_rover) * Frr_simple
    
#     return Frr

# @nb.njit()
# def __F_rolling(omega, terrain_angle, Crr, m, g, r, Ng):
#     # initialize
#     Frr = np.zeros_like(omega)

#     for ii in range(len(omega)):
#         v_rover = r*omega[ii]/Ng
#         Fn = m*g*math.cos(math.radians(terrain_angle[ii])) # normal force
#         Frr_simple = -Crr*Fn # simple rolling resistance
#         Frr[ii] = math.erf(40*v_rover) * Frr_simple
    
#     return Frr

# def F_rolling(omega, terrain_angle, rover, planet, Crr):
#     """
#     Inputs:           omega:  numpy array     Motor shaft speed [rad/s]
#               terrain_angle:  numpy array     Array of terrain angles [deg]
#                       rover:  dict            Data structure specifying rover 
#                                               parameters
#                     planet:  dict            Data dictionary specifying planetary 
#                                               parameters
#                         Crr:  scalar          Value of rolling resistance coefficient
#                                               [-]
    
#     Outputs:           Frr:  numpy array     Array of forces [N]
#     """
    
#     # Check that the first input is a scalar or a vector
#     if (type(omega) != int) and (type(omega) != float) and (not isinstance(omega, np.ndarray)):
#         raise Exception('First input must be a scalar or a vector. If input is a vector, it should be defined as a numpy array.')
#     elif not isinstance(omega, np.ndarray):
#         omega = np.array([omega],dtype=float) # make the scalar a numpy array
#     elif len(np.shape(omega)) != 1:
#         raise Exception('First input must be a scalar or a vector. Matrices are not allowed.')
        
#     # Check that the second input is a scalar or a vector
#     if (type(terrain_angle) != int) and (type(terrain_angle) != float) and (not isinstance(terrain_angle, np.ndarray)):
#         raise Exception('Second input must be a scalar or a vector. If input is a vector, it should be defined as a numpy array.')
#     elif not isinstance(terrain_angle, np.ndarray):
#         terrain_angle = np.array([terrain_angle],dtype=float) # make the scalar a numpy array
#     elif len(np.shape(terrain_angle)) != 1:
#         raise Exception('Second input must be a scalar or a vector. Matrices are not allowed.')
        
#     # Check that the first two inputs are of the same size
#     if len(omega) != len(terrain_angle):
#         raise Exception('First two inputs must be the same size')
    
#     # Check that values of the second input are within the feasible range  
#     if max([abs(x) for x in terrain_angle]) > 75:    
#         raise Exception('All elements of the second input must be between -75 degrees and +75 degrees')
        
#     # Check that the third input is a dict
#     if type(rover) != dict:
#         raise Exception('Third input must be a dict')
        
#     # Check that the fourth input is a dict
#     if type(planet) != dict:
#         raise Exception('Fourth input must be a dict')
        
#     # Check that the fifth input is a scalar and positive
#     if (type(Crr) != int) and (type(Crr) != float):
#         raise Exception('Fifth input must be a scalar')
#     if Crr <= 0:
#         raise Exception('Fifth input must be a positive number')
        
#     # Main Code
#     m = get_mass(rover)
#     g = planet['g']
#     r = rover['wheel_assembly']['wheel']['radius']
#     Ng = get_gear_ratio(rover['wheel_assembly']['speed_reducer'])
    
#     # v_rover = r*omega/Ng
    
#     # # Fn = np.array([m*g*math.cos(math.radians(x)) for x in terrain_angle],dtype=float) # normal force
#     # Fn = m*g*np.cos(np.radians(terrain_angle)) # normal force
#     # Frr_simple = -Crr*Fn # simple rolling resistance
    
#     # Frr = np.array([math.erf(40*v_rover[ii]) * Frr_simple[ii] for ii in range(len(v_rover))], dtype = float)
#     Frr = __F_rolling(omega, terrain_angle, Crr, m, g, r, Ng)
    
#     return Frr


# def motorW(v, rover):
#     '''
#     Compute the rotational speed of the motor shaft [rad/s] given the translational velocity of the rover and the rover
#     dictionary.
#     This function should be “vectorized” such that if given a vector of rover velocities it returns a vector the same size
#     containing the corresponding motor speeds.

#     Inputs
#     ------
#     v: scalar or 1D numpy array
#         Rover translational velocity [m/s]
#     rover: dict
#         Data structure containing rover parameters

#     Outputs
#     -------
#     w: scalar or 1D numpy array
#         Motor speed [rad/s]
#         Return argument should match type/size of input
#     '''

#     # check that the velocity argument is a scalar or numpy array
#     if not isinstance(v, (np.ndarray, numbers.Number)):
#         raise TypeError('1st argument (v) must be a scalar or a numpy array')
#     if isinstance(v, np.ndarray) and v.ndim != 1:
#         raise TypeError(f'1st argument (v) must be a scalar or a vector, not {type(v)} Matricies are not allowed. v.ndim == {np.shape(v)}')

#     try:
#         gear_ratio = get_gear_ratio(rover['wheel_assembly']['speed_reducer'])  # ratio of gearbox input shaft speed to output shaft speed
#         wheel_w = v / rover['wheel_assembly']['wheel']['radius']  # rotational velocity of the wheel [rad/s]
#     except KeyError as e:
#         raise KeyError(f'Invalid rover dictionary, could not find key: {e}')

#     motor_w = wheel_w * gear_ratio

#     return motor_w

# def quick_motorW_factory(rover):
#     gear_ratio = get_gear_ratio(rover['wheel_assembly']['speed_reducer'])  # ratio of gearbox input shaft speed to output shaft speed
#     wheel_radius = rover['wheel_assembly']['wheel']['radius']
#     def quick_motorW(v):
#         return v / wheel_radius * gear_ratio
#     return quick_motorW

# def __rover_dynamics():
#     """
#     Inputs:           None
    
#     Outputs:          None
#     """
    
#     # Main Code
#     pass


def rover_dynamics(t, y, rover, planet, experiment, alpha_fun = None):
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
    if not isinstance(t, numbers.Number):
        raise Exception (f"1st argument (t) must be a scalar, not {type(t)}")
    # check that y is a np array
    if not isinstance(y, np.ndarray):
        raise Exception(f"2nd argument (y) must be a numpy array, not {type(y)}")
    # check that the vector is 1D
    if y.ndim != 1:
        raise TypeError('2nd argument (y) must be a vector. Matricies are not allowed.')
    #check that y only has two elements
    if len(y) != 2:
        raise Exception("2nd argument (y) must only contain two elements, acceleration and velocity")

    # check that rover, planet, and experiment are dicitonaries 
    if (type(rover) != dict) or (type(planet) != dict) or (type(experiment) != dict):
        raise Exception("Third, fourth and fifth arguments must be dictionaries")
        
    # define alpha_fun
    if alpha_fun is None:
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

# def quick_dynamics_factory(rover, planet, experiment):
    # '''
    # returns a function which takes in a time and state vector and returns the derivative of the state vector
    # '''
    # alpha_fun = interp1d(experiment['alpha_dist'], experiment['alpha_deg'], kind = 'cubic', fill_value = 'extrapolate')
    # mass = get_mass(rover)
    # quick_motorW = quick_motorW_factory(rover)
    # def quick_dynamics(t, y):
    #     terrain_angle = float(alpha_fun(y[1]))
    #     velocity = float(y[0])
    #     a = np.float64((1/mass) * F_net(quick_motorW(velocity), terrain_angle, rover, planet, experiment['Crr']))
    #     dydt = np.array([a, y[0]])
    #     return dydt
    # return quick_dynamics


def simulate_rover(rover: dict, planet: dict, experiment: dict, end_event: dict=None, quick_rover_dynamics = None, expand: bool=False) -> dict:
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
    # check if rover is dict
    if not isinstance(rover, dict):
        raise TypeError(f"1st argument (rover) must be a dictionary, not {type(rover)}")
    if not isinstance(planet, dict):
        raise TypeError(f"2nd argument (planet) must be a dictionary, not {type(planet)}")
    if not isinstance(experiment, dict):
        raise TypeError(f"3rd argument (experiment) must be a dictionary, not {type(experiment)}")
    if (end_event is not None) and (not isinstance(end_event, dict)):
        raise TypeError(f"optional kwarg end_event must be a dictionary, not {type(end_event)}")
    if not isinstance(expand, bool):
        raise TypeError(f"optional kwarg expand must be a boolean, not {type(expand)}")

    # allow end_event to override experiment['end_event'] if it is provided
    if end_event is None:
        end_event = experiment['end_event']

    time_span = experiment['time_range']
    y0 = experiment['initial_conditions']

    # functools.partial allows us to pass in arguments to the function that will be called by the integrator
    # you can also use lambda functions, but this is more readable to me
    alpha_fun = interp1d(experiment['alpha_dist'], experiment['alpha_deg'], kind = 'cubic', fill_value = 'extrapolate')
    rover_dynamics_partial = ft.partial(rover_dynamics, rover=rover, planet=planet, experiment=experiment, alpha_fun=alpha_fun)
    sol = solve_ivp(rover_dynamics_partial, time_span, y0, method='BDF', events=end_of_mission_event(end_event), dense_output=expand, rtol=1e-10, atol=1e-10)
    t = sol.t
    y = sol.y

    # calculate and store telemetry data
    telemetry = {
        'time': t,
        'completion_time': t[-1],
        'velocity': y[0],
        'position': y[1],
        'distance_traveled': y[1][-1] - y[1][0],  # final position - initial position
        'max_velocity': np.max(y[0]),
        'average_velocity': (y[1][-1] - y[1][0])/t[-1],  # distance traveled / time elapsed
        'power': mechpower(y[0], rover),
        'battery_energy': battenergy(t, y[0], rover),
        'energy_per_distance': battenergy(t, y[0], rover)/(y[1][-1] - y[1][0])  # total energy / distance traveled
    }
    # allow the function caller to access the solution object if desired
    if expand:
        telemetry['sol'] = sol
    # make a copy of the rover dictionary so we don't modify the original
    rover = rover.copy()
    rover['telemetry'] = telemetry
    return rover