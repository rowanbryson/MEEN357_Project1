
# This is a template for the Marvin dictionary
# I just converted the text from the project description into a dictionary
rover = {
    'wheel_assembly': {
        'wheel': {
            'radius': 0.30,  # Radius of drive wheel [m]
            'mass': 1.0,  # Mass of one drive wheel [kg]
        },
        'speed_reducer': {
            'type': "reverted",  # String of text defining the type of speed reducer. For Project Phase 1, the only valid entry is “reverted”.
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
planet = {
    'g': 3.72  # Acceleration due to gravity [m/s^2]
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

'''
Yo -Andrew
'''
