# this is replacing the old "MARVIN_DICT"
rover_1_dict = {
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