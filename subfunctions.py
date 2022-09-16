
# This is a template for the Marvin dictionary
# I just converted the text from the project description into a dictionary
MARVIN_DICT_TEMPLATE = {
    'rover': {
        'wheel_assembly': {
            'wheel': {
                'radius': None,  # Radius of drive wheel [m]
                'mass': None,  # Mass of one drive wheel [kg]
            },
            'speed_reducer': {
                'type': None,  # String of text defining the type of speed reducer. For Project Phase 1, the only valid entry is “reverted”.
                'diam_pinion': None,  # Diameter of pinion [m]
                'diam_gear': None,  # Diameter of gear [m]
                'mass': None,  # Mass of speed reducer assembly [kg]
            },
            'motor': {
                'torque_stall': None,  # Motor stall torque [Nm]
                'torque_noload': None,  # Motor no-load torque [Nm]
                'speed_noload': None,  # Motor no-load speed [rad/s]
                'mass': None,  # Motor mass [kg]
            }
        },
        'chassis': {
            'mass': None  # Mass of chassis [kg]
        },
        'science_payload': {
            'mass': None  # Mass of science payload [kg]
        },
        'power_subsys': {
            'mass': None  # Mass of power subsystem [kg]
        },
    },
    'planet': {
        'g': None  # Acceleration due to gravity [m/s^2]
    }
}


def get_mass():
    pass

def get_gear_ratio():
    pass

def tau_dcmotor():
    pass

def F_drive():
    pass

def F_rolling():
    pass

def F_net():
    pass