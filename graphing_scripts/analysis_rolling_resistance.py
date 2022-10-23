# allow imports from parent directory
import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from subfunctions import *
import define_rovers
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import root_scalar
from scipy.optimize import root

# this is a band aid to make the code work, the new subfunctions.py file doesn't define MARVIN_DICT
# we might eventually want to change this to a function that takes in a rover name
MARVIN_DICT = {'rover': define_rovers.rover1(), 'planet': define_rovers.planet1()}

def main(save_plots=True):
    
    #Initialize Variables
    terrain_angle = 0
    Crr_array = np.linspace(0.01, 0.4, num=25)
    v_max = np.zeros(len(Crr_array), dtype=float)
    omega_nl = MARVIN_DICT['rover']['wheel_assembly']['motor']['speed_noload']
    wheel_radius = MARVIN_DICT['rover']['wheel_assembly']['wheel']['radius']
    rover = MARVIN_DICT['rover']
    planet = MARVIN_DICT['planet']
    gear_ratio = get_gear_ratio(rover['wheel_assembly']['speed_reducer'])

    for ii in range(len(Crr_array)):
        fun = lambda omega: F_net(np.array([omega]), np.array([float(terrain_angle)]), rover, planet, float(Crr_array[ii]))
        sol = root_scalar(fun, method='bisect', bracket=[0, omega_nl])
        omega_max = sol.root
        v_max[ii] = omega_max / gear_ratio * rover['wheel_assembly']['wheel']['radius']

    #Plot acceleration vs. Crr
    plt.plot(Crr_array, v_max)
    plt.title('Max Velocity vs. Rolling Resistance Coefficient')
    plt.ylabel('V Max (m/s)')
    plt.xlabel('Crr')
    if save_plots:
        plt.savefig('plots/analysis_rolling_resistance.png')
    plt.show()


if __name__ == '__main__':
        main()