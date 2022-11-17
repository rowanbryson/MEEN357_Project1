"""###########################################################################
#   This file initializes a rover structure for testing/grading
#
#   Created by: Jonathan Weaver-Rosen
#   Last Modified: 25 June 2021
###########################################################################"""

import numpy as np

def define_rover_1():
    # Initialize Rover dict for testing
    wheel = {'radius':0.30,
             'mass':1}
    speed_reducer = {'type':'reverted',
                     'diam_pinion':0.04,
                     'diam_gear':0.07,
                     'mass':1.5}
    motor = {'torque_stall':170,
             'torque_noload':0,
             'speed_noload':3.80,
             'mass':5.0}
    
    
    # phase 2 add ##############################
    motor['effcy_tau'] = np.array([0, 10, 20, 40, 75, 165])
    motor['effcy']     = np.array([0,.60,.75,.73,.55, .05])
    #############################################
    
    
    chassis = {'mass':659}
    science_payload = {'mass':75}
    power_subsys = {'mass':90}
    
    wheel_assembly = {'wheel':wheel,
                      'speed_reducer':speed_reducer,
                      'motor':motor}
    
    rover = {'wheel_assembly':wheel_assembly,
             'chassis':chassis,
             'science_payload':science_payload,
             'power_subsys':power_subsys}
    
    planet = {'g':3.72}
    
    # end_event = {}
    # end_event['max_distance'] = 50
    # end_event['max_time'] = 2000
    # end_event['min_velocity'] = 0.01
    # end_event['on_surface'] = 'False'
    
    # experiment = {}
    # experiment['crr'] = 0.1
    # experiment['time_range'] = np.array([0, 2000])
    # experiment['initial_conditions'] = np.array([[0.2],[0]])
    # experiment['alpha_deg']  = np.array([0, -3, 3, 10, 5])
    # experiment['alpha_dist'] = np.array([0, 5, 15, 30, 50])
    
    # return everything we need
    return rover, planet

def define_rover_2():
    # Initialize Rover dict for testing
    wheel = {'radius':0.30,
             'mass':2} 
    speed_reducer = {'type':'reverted',
                     'diam_pinion':0.04,
                     'diam_gear':0.06,
                     'mass':1.5}
    motor = {'torque_stall':180,
             'torque_noload':0,
             'speed_noload':3.70,
             'mass':5.0}
    
    
    # # phase 2 add ##############################
    # motor['effcy_tau'] = [0, 10, 20, 40, 75, 165]
    # motor['effcy']     = [0,.60,.75,.73,.55, .05]
    # #############################################
    
    
    chassis = {'mass':659}
    science_payload = {'mass':75}
    power_subsys = {'mass':90}
    
    wheel_assembly = {'wheel':wheel,
                      'speed_reducer':speed_reducer,
                      'motor':motor}
    
    rover = {'wheel_assembly':wheel_assembly,
             'chassis':chassis,
             'science_payload':science_payload,
             'power_subsys':power_subsys}
    
    planet = {'g':3.72}
    
    # end_event = {}
    # end_event['max_distance'] = 50
    # end_event['max_time'] = 2000
    # end_event['min_velocity'] = 0.01
    # end_event['on_surface'] = 'False'
    
    # experiment = {}
    # experiment['crr'] = 0.1
    # experiment['time_range'] = [0, 2000]
    # experiment['initial_conditions'] = [[0.2][0]]
    # experiment['alpha_deg']  = [0, -3, 3, 10, 5]
    # experiment['alpha_dist'] = [0, 5, 15, 30, 50]
    
    # return everything we need
    return rover, planet

def define_rover_3():
    # Initialize Rover dict for testing
    wheel = {'radius':0.30,
             'mass':2} 
    speed_reducer = {'type':'standard',
                     'diam_pinion':0.04,
                     'diam_gear':0.06,
                     'mass':1.5}
    motor = {'torque_stall':180,
             'torque_noload':0,
             'speed_noload':3.70,
             'mass':5.0}
    
    
    # # phase 2 add ##############################
    # motor['effcy_tau'] = [0, 10, 20, 40, 75, 165]
    # motor['effcy']     = [0,.60,.75,.73,.55, .05]
    # #############################################
    
    
    chassis = {'mass':659}
    science_payload = {'mass':75}
    power_subsys = {'mass':90}
    
    wheel_assembly = {'wheel':wheel,
                      'speed_reducer':speed_reducer,
                      'motor':motor}
    
    rover = {'wheel_assembly':wheel_assembly,
             'chassis':chassis,
             'science_payload':science_payload,
             'power_subsys':power_subsys}
    
    planet = {'g':3.72}
    
    # return everything we need
    return rover, planet


def define_rover_4():
    # Initialize Rover dict for testing
    wheel = {'radius':0.20,
             'mass':2} 
    speed_reducer = {'type':'reverted',
                     'diam_pinion':0.04,
                     'diam_gear':0.06,
                     'mass':1.5}
    motor = {'torque_stall':165,
             'torque_noload':0,
             'speed_noload':3.85,
             'mass':5.0}
    
    # phase 2 add ##############################
    motor['effcy_tau'] = np.array([0, 10, 20, 40, 75, 170])
    motor['effcy']     = np.array([0,.60,.75,.73,.55, .05])
    #############################################
    
    
    chassis = {'mass':674}
    science_payload = {'mass':80}
    power_subsys = {'mass':100}
    
    wheel_assembly = {'wheel':wheel,
                      'speed_reducer':speed_reducer,
                      'motor':motor}
    
    rover = {'wheel_assembly':wheel_assembly,
             'chassis':chassis,
             'science_payload':science_payload,
             'power_subsys':power_subsys}
    
    planet = {'g':3.72}
    
    end_event = {'max_distance' : 10,
                 'max_time' : 10000,
                 'min_velocity' : 0.01}
    
    # return everything we need
    return rover, planet #, end_event