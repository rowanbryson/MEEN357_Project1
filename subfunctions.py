import numpy as np

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
    if not isinstance(rover, dict):
        raise TypeError('Input argument rover must be a dictionary.')

    # Wheel Assembly Calculation
    wheel = rover['wheel_assembly']['wheel']['mass']
    speed_reducer = rover['wheel_assembly']['speed_reducer']['mass']
    motor = rover['wheel_assembly']['motor']['mass']
    wheel_assembly = wheel + speed_reducer + motor

    # Other Masses
    chassis = rover['chassis']['mass']
    science_payload = rover['science_payload']['mass']
    power_subsys = rover['power_subsys']['mass']

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

    #  TODO add input checking

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
    Fd = tau * gear_ratio / wheel_radius
    return Fd


def F_gravity(terrain_angle: np.ndarray, rover: dict, planet: dict):
    pass


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
    pass


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

    # Force Calculations
    drive = F_drive(omega, rover)
    rolling = F_rolling(omega, terrain_angle, rover, planet, Crr)
    gravity = F_gravity(terrain_angle, rover, planet)

    # Force in the direction of motion
    return drive - rolling + gravity*np.sin(terrain_angle)



'''
Yo -Andrew
'''
