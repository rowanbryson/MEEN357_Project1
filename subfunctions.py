import numpy as np
# This is a template for the Marvin dictionary
# I just converted the text from the project description into a dictionary
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
    pass


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
    pass


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
    pass


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
    pass


'''
Yo -Andrew
'''
