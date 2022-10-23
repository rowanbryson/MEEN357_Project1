from subfunctions import *
import matplotlib.pyplot as plt
from scipy.optimize import root_scalar
from scipy.optimize import root
import numpy as np

def main(save_plots=True):
    Crr = .2
    slope_list_deg = np.linspace(-10, 35, 25)
    v_max = np.zeros(len(slope_list_deg), dtype=float)
    omega_nl = MARVIN_DICT['rover']['wheel_assembly']['motor']['speed_noload']
    wheel_radius = MARVIN_DICT['rover']['wheel_assembly']['wheel']['radius']
    rover = MARVIN_DICT['rover']
    planet = MARVIN_DICT['planet']
    gear_ratio = get_gear_ratio(rover['wheel_assembly']['speed_reducer'])

    for ii in range(len(slope_list_deg)):
        fun = lambda omega: F_net(np.array([omega]), np.array([float(slope_list_deg[ii])]), rover, planet, Crr)
        sol = root_scalar(fun, method='bisect', bracket=[0, omega_nl])
        omega_max = sol.root
        v_max[ii] = omega_max / gear_ratio * rover['wheel_assembly']['wheel']['radius']

    plt.plot(slope_list_deg, v_max)
    plt.title("Slope [deg] vs Maximum Velocity [m/s]")
    plt.xlabel("Slope [deg]")
    plt.ylabel("Maximum Velocity [m/s]")
    if save_plots:
        plt.savefig('plots/analysis_terrain_slope.png')
    plt.show()

if __name__ == '__main__':
        main()