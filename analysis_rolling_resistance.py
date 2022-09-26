from subfunctions import *
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import root_scalar
from scipy.optimize import root

def main():
    
    #Initialize Variables
    terrain_angle = 0
    Crr_array = np.linspace(0.01, 0.4, num=25)
    v_max = np.zeros(len(Crr_array), dtype=float)
    omega_nl = MARVIN_DICT['rover']['wheel_assembly']['motor']['speed_noload']
    wheel_radius = MARVIN_DICT['rover']['wheel_assembly']['wheel']['radius']
    rover = MARVIN_DICT['rover']
    planet = MARVIN_DICT['planet']

    for ii in range(len(Crr_array)):
        fun = lambda omega: F_net(np.array([omega]), np.array([float(terrain_angle)]), rover, planet, Crr_array[ii])
        sol = root_scalar(fun, method='bisect', bracket=[0, omega_nl])
        omega_max = sol.root
        v_max[ii] = wheel_radius * 6 * omega_max

    #Plot acceleration vs. Crr
    plt.plot(Crr_array, v_max)
    plt.title('Max Velocity vs. Rolling Resistance Coefficient')
    plt.ylabel('V Max (m/s)')
    plt.xlabel('Crr')
    plt.show()


if __name__ == '__main__':
        main()