import numpy as np
import math

# The following dict is the MARVIN_DICT specified in the project description
# not sure if it will live here forever but it's here for now
MARVIN_DICT = {
    'rover': {
        'wheel_assembly': {
            'wheel': {
                'radius': 0.30,  # Radius of drive wheel [m]
                'mass': 1.0,  # Mass of one drive wheel [kg]
            },
            'speed_reducer': {
                # String of text defining the type of speed reducer. For Project Phase 1, the only valid entry is “reverted”.
                'type': 'reverted',
                'diam_pinion': 0.04,  # Diameter of pinion [m]
                'diam_gear': 0.07,  # Diameter of gear [m]
                'mass': 1.5,  # Mass of speed reducer assembly [kg]
            },
            'motor': {
                'torque_stall': 170,  # Motor stall torque [Nm]
                'torque_noload': 0,  # Motor no-load torque [Nm]
                'speed_noload': 3.80,  # Motor no-load speed [rad/s]
                'mass': 5.0,  # Motor mass [kg]
            }
        },
        'chassis': {
            'mass': 659  # Mass of chassis [kg]
        },
        'science_payload': {
            'mass': 75  # Mass of science payload [kg]
        },
        'power_subsys': {
            'mass': 90  # Mass of power subsystem [kg]
        },
    },
    'planet': {
        'g': 3.72  # Acceleration due to gravity [m/s^2]
    }
}


def get_mass(rover: dict) -> float:
    '''
    This function computes rover mass in kilograms. It accounts for the chassis, power subsystem, science payload,
    and six wheel assemblies, which itself is comprised of a motor, speed reducer, and the wheel itself.

    Parameters
    ----------
    rover : dict
        Data structure containing rover parameters

    Returns
    -------
    mass : scalar
        Rover mass in [kg]

    '''
    # Input Validation

    # (a) input is dict
    if not isinstance(rover, dict):
        raise Exception('Input argument must be a dictionary.')

    # Wheel Assembly Calculation
    wheel = rover['wheel_assembly']['wheel']['mass']
    speed_reducer = rover['wheel_assembly']['speed_reducer']['mass']
    motor = rover['wheel_assembly']['motor']['mass']
    wheel_assembly = wheel + speed_reducer + motor

    # Other Masses
    chassis = rover['chassis']['mass']
    science_payload = rover['science_payload']['mass']
    power_subsys = rover['power_subsys']['mass']

    # Calculation / Output
    return 6 * wheel_assembly + chassis + science_payload + power_subsys


def get_gear_ratio(speed_reducer: dict) -> float:
    '''
    This function computes the gear ratio of the speed reducer.
    In later project phases, you will extend this to work for various types of speed reducers. For now, it needs to work
    only with the simple reverted gear set described in Section 2.2.

    Parameters
    ----------
    speed_reducer : dict
        Data structure specifying speed reducer parameters
    
    Returns
    -------
    gear_ratio : scalar
        Speed ratio from input pinion shaft to output gear shaft. Unitless.
    '''

    # Input Validation

    # (a) input is dict
    if not isinstance(speed_reducer, dict):
        raise Exception('Input argument must be a dictionary.')

    # Comparison
    if not speed_reducer['type'].lower() == 'reverted':
        raise Exception('get_gear_ratio does not currently support speed reducer types other than reverted')

    # Calculation / Output
    return (speed_reducer['diam_gear'] / speed_reducer['diam_pinion'])**2


def tau_dcmotor(omega: np.ndarray, motor: dict) -> np.ndarray:
    '''
    This function returns the motor shaft torque in Nm given the shaft speed in rad/s and the motor specifications
    structure (which defines the no-load speed, no-load torque, and stall speed, among other things.

    This function should operate in a “vectorized” manner, meaning that if given a vector of motor shaft speeds, it
    returns a vector of the same size consisting of the corresponding motor shaft torques.

    Parameters
    ----------
    omega: numpy array
        Motor shaft speed [rad/s]
    
    motor: dict
        Data structure specifying motor parameters

    Returns
    -------
    tau: numpy array
        Torque at motor shaft [Nm]. Return argument is same size as first input argument.
    '''

    #### INPUT CHECKING ####

    # if input is arraylike but not a numpy array, convert to a numpy array
    if not isinstance(omega, np.ndarray):
        try:
            if isinstance(omega, (int, float)):
                omega = np.array([omega], dtype=np.float64)
            else:
                omega = np.array(omega, dtype=np.float64)
        except:
            raise TypeError('Input argument omega must be arraylike and convertible to a numpy array of floats.')
    # check that motor dict is valid format
    if not isinstance(motor, dict):
        raise TypeError('Input argument motor must be a dictionary.')
    if not all([key in motor for key in ['torque_stall', 'torque_noload', 'speed_noload']]):
        raise ValueError('Input argument motor must contain keys "torque_stall", "torque_noload", and "speed_noload".')
    if not all([isinstance(motor[key], (int, float)) for key in ['torque_stall', 'torque_noload', 'speed_noload']]):
        raise TypeError('"torque_stall", "torque_noload", and "speed_noload" must be numeric scalars.')

    #### FUNCTION BODY ####

    # initialize tau as a numpy array of the same size as omega with float values of 0
    tau = np.zeros_like(omega, dtype=float)
    # define the function to calculate the torque under normal conditions
    main_tau_func = lambda omega: motor['torque_stall'] - (motor['torque_stall'] - motor['torque_noload']) * (omega / motor['speed_noload'])

    # iterate through omega and calculate the torque at each speed
    for i in range(len(omega)):
        # see the pdf for the logic behind this
        # 'omega > noload_speed' and 'omega < 0' are special cases
        if 0 <= float(omega[i]) <= motor['speed_noload']:
            tau[i] = main_tau_func(omega[i])
        elif omega[i] > motor['speed_noload']:
            tau[i] = 0
        else:  # omega < 0
            tau[i] = motor['torque_stall']
    return tau


def F_drive(omega: np.ndarray, rover: dict) -> np.ndarray:
    '''
    This function computes the combined drive force, in Newtons, acting on the rover due to all six wheels. This
    force is a function of the motor shaft speed, and the properties of the motor, speed reducer, and drive track (all
    defined in the rover dict).
    This function should be “vectorized” such that if given a vector of motor shaft speeds, it returns a vector of the
    same size consisting of the corresponding forces.

    Parameters
    ----------
    omega: numpy array
        Array of motor shaft speeds [rad/s]
    
    rover: dict
        Data structure specifying rover parameters

    Returns
    -------
    Fd: numpy array
        Array of drive forces [N]
    '''
    gear_ratio = get_gear_ratio(rover['wheel_assembly']['speed_reducer'])
    wheel_radius = rover['wheel_assembly']['wheel']['radius']

    tau = tau_dcmotor(omega, rover['wheel_assembly']['motor'])
    # multiply by 6 because there are 6 wheels
    Fd = 6 * tau * gear_ratio / wheel_radius

    # Fd = np.zeros_like(omega, dtype=float)
    # for i in range(len(omega)):
    #     Fd[i] = 6 * tau[i] * gear_ratio / wheel_radius

    return Fd


def F_gravity(terrain_angle: np.ndarray, rover: dict, planet: dict):
    '''
    This function computes the component of force due to gravity, in Newtons, acting in the direction of rover
    translation. This force is a function of the angle the terrain makes with the horizon (in degrees) and the total mass
    of the rover (in kg).
    This function should be “vectorized” such that if given a vector of terrain angles, it returns a vector of the same
    size consisting of the corresponding forces.

    Parameters
    ----------
    terrain_angle: numpy array
        Array of terrain angles [deg]
            !! in degrees !!

    rover: dict
        Data structure containing rover parameters

    planet: dict
        Data structure containing gravity parameter

    Returns
    -------
    Fgt: numpy array
        Array of gravity forces [N]
    '''

    if not isinstance(rover, dict):
        raise Exception("rover must be a dictionary")
    if not isinstance(planet, dict):
        raise Exception("planet must be a dictionary")

    if isinstance(terrain_angle, (int, float)):
        terrain_angle = np.array([terrain_angle], dtype=np.float64)
    elif isinstance(terrain_angle, list):
        terrain_angle = np.array(terrain_angle, dtype=np.float64)
    if not isinstance(terrain_angle, np.ndarray):
        raise Exception("terrain_angle must be an arraylike or a scalar")
    
    # check that all values in terrain_angle are between -75 and 75
    if not np.all(np.logical_and(terrain_angle >= -75, terrain_angle <= 75)):
        raise Exception("all values in terrain_angle must be between -75 and 75")

    #adds the gravity force to empty array for each terrain angle 
    Fgt = np.zeros_like(terrain_angle)
    for i in range(len(terrain_angle)):

        terrain_angle_rad = np.deg2rad(terrain_angle[i])
        #F_gravity = (-3.72)*cos(terrainangle)*(total mass of rover)
        gravity_force = (planet['g'])*np.sin(terrain_angle_rad)*get_mass(rover)*(-1)

        Fgt[i] = gravity_force
    return Fgt


def F_rolling(omega: np.ndarray, terrain_angle: np.ndarray, rover: dict, planet: dict, Crr: float) -> np.ndarray:
    '''
    This function computes the component of force due to rolling resistance, in Newtons, acting in the direction of
    rover translation. This force is a function of the angle the terrain makes with the horizon (in degrees), the total
    mass of the rover (in kg), and the rolling resistance coefficient (which is a unitless constant that depends on
    properties of both the wheels and the ground).
    This function should be “vectorized” such that if given a vector of terrain angles, it returns a vector of the same
    size consisting of the corresponding forces.
    This function computes the rolling resistance summed over all six wheels. Assume that 1/6th the rover normal force
    acts on each wheel.

    Parameters
    ----------
    omega: numpy array
        Array of motor shaft speeds [rad/s]

    terrain_angle: numpy array
        Array of terrain angles [deg]

    rover: dict
        Data structure containing rover parameters

    planet: dict
        Data structure containing planet gravity parameter

    Crr: scalar
        Value of rolling resistance coefficient [-]

    Returns
    -------
    Frr: numpy array
        Array of rolling resistance forces [N]
    '''
    #Check for inital conditions

    if not isinstance(rover, dict) or not isinstance(planet, dict):
        raise Exception("The rover and planet attributes must be dictionaries.")
    if not isinstance(terrain_angle, np.ndarray) or not isinstance(omega, np.ndarray):
        raise Exception("The terrain angle must be in an numpy array.")
    if not np.issubdtype(omega.dtype, np.number):
        raise Exception("The omega array must be of type float64.")
    if not np.issubdtype(terrain_angle.dtype, np.number):
        raise Exception(f"The terrain angle array must be of type float64 or int64, not {terrain_angle.dtype}.")
    if not omega.shape == terrain_angle.shape:
        raise Exception("The omega and terrain angle arrays must be the same size.")
    if np.any(terrain_angle < -75) or np.any(terrain_angle > 75):
        raise Exception("All terrain angles must be between -75 and 75 degrees.")
    if (not isinstance(Crr, int) and not isinstance(Crr, float)) or (Crr < 0):
        raise Exception("The value of Crr must be a positive scalar value.")
    
    Frr = np.zeros_like(terrain_angle, dtype=float)

    gear_ratio = get_gear_ratio(rover['wheel_assembly']['speed_reducer'])
    wheel_radius = rover['wheel_assembly']['wheel']['radius']

    for i in range(len(terrain_angle)):
        wheel_velocity = (omega[i] / gear_ratio) * (wheel_radius)
        normal_force = get_mass(rover) * planet['g'] * np.cos(np.deg2rad(terrain_angle[i]))
        force_rolling = -Crr * normal_force * math.erf(40 * wheel_velocity)
        Frr[i] = force_rolling

    return Frr


def F_net(omega: np.ndarray, terrain_angle: np.ndarray, rover: dict, planet: dict, Crr: float) -> np.ndarray:
    '''
    This function computes the total force, in Newtons, acting on the rover in the direction of its motion. This is a
    function of the motor shaft speed, terrain angle, rover properties, and the coefficient of rolling resistance.

    This function should be “vectorized” such that if given same-sized vectors of motor shaft speeds and terrain angles,
    it returns a vector of the same size consisting of the corresponding forces.

    Parameters
    ----------
    omega: numpy array
        Array of motor shaft speeds [rad/s]

    terrain_angle: numpy array
        Array of terrain angles [deg]

    rover: dict
        Data structure containing rover parameters

    planet: dict
        Data structure containing planet gravity parameter

    Crr: scalar
        Value of rolling resistance coefficient [-]

    Returns
    -------
    Fnet: numpy array
        Array of net forces [N]
    '''

    # Input Validation

    # (a) first and second arguments are a scalar or vector
    if isinstance(omega, (int, float)):
        omega = np.array([omega])
    elif isinstance(omega, list):
        omega = np.array(omega)
    if not isinstance(omega, np.ndarray):
        raise Exception('omega should be array or scalar.')

    if isinstance(terrain_angle, (int, float)):
        terrain_angle = np.array([terrain_angle])
    elif isinstance(terrain_angle, list):
        terrain_angle = np.array(terrain_angle)
    if not isinstance(terrain_angle, np.ndarray):
        raise Exception('terrain_angle should be array or scalar.')

    # (a) first and second arguments are equal size
    if not (omega.size == terrain_angle.size):
        raise Exception('omega and terrain_angle are not same size')

    # (b) elements of second are between -75 and 75
    outside_range = False
    for i in terrain_angle:
        if i < -75 or i > 75:
            outside_range = True
    if outside_range is True:
        raise Exception("terrain_angle must be between -75 and 75 degrees")

    # (c) third and fourth are dict
    if not isinstance(rover, dict):
        raise Exception("rover should be a dict")
    if not isinstance(planet, dict):
        raise Exception("planet should be a dict")

    # (d) fifth is positive scalar
    if not isinstance(Crr, (int, float)):
        raise Exception("Crr should be a scalar")
    if not Crr > 0:
        raise Exception("Crr should be positive")


    # Force Calculations
    drive = F_drive(omega, rover)
    rolling = F_rolling(omega, terrain_angle, rover, planet, Crr)
    gravity = F_gravity(terrain_angle, rover, planet)

    return drive + rolling + gravity
