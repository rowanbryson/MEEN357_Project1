import sys
sys.path.append('../')
from subfunctions import *
import define_rovers
import matplotlib.pyplot as plt
from scipy.optimize import root_scalar
import numpy as np

# the code works great, don't worry about this
import warnings
warnings.filterwarnings("ignore")

# this is a band aid to make the code work, the new subfunctions.py file doesn't define MARVIN_DICT
# we might eventually want to change this to a function that takes in a rover name
MARVIN_DICT = define_rovers.default_data_dict

def main(save_plots=False):
    # initialize domain of interest
    Crr_array = np.linspace(0.01, 0.4, 25)
    slope_array_deg = np.linspace(-10, 35, 25)

    # initialize empty arrays to store results
    CRR, SLOPE = np.meshgrid(Crr_array, slope_array_deg)
    VMAX = np.zeros(CRR.shape, dtype=float)

    # initialize rover and planet
    rover = MARVIN_DICT['rover']
    planet = MARVIN_DICT['planet']
    gear_ratio = get_gear_ratio(rover['wheel_assembly']['speed_reducer'])

    # choose initial guess for omega_max
    no_load_omega = MARVIN_DICT['rover']['wheel_assembly']['motor']['speed_noload']
    x0_guess = no_load_omega * 0.9
    x1_guess = no_load_omega * 0.1

    # calculate maximum velocity
    for i in range(len(Crr_array)):
        for j in range(len(slope_array_deg)):
            # this defines a function of Net force in terms of omega
            # the root of this function is omega_max
            # which is the maximum motor shaft speed
            F_net_from_omega = lambda omega: F_net(omega, float(slope_array_deg[j]), rover, planet, Crr_array[i])

            # I chose secant because it doesn't require a bracket
            # using 0 to noload_speed doesn't work because the function doesn't seem to be monotonic
            sol = root_scalar(F_net_from_omega, method='secant', x0=x0_guess, x1=x1_guess)
            omega_max = sol.root
            VMAX[i,j] = omega_max / gear_ratio * rover['wheel_assembly']['wheel']['radius']
            

    # plot in 3D
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(CRR, SLOPE, VMAX, cmap='viridis', edgecolor='none')
    ax.set_title("Max Velocity vs Coeff of Rolling Resistance (Crr) and Terrain Slope")
    ax.set_xlabel("Crr")
    ax.set_ylabel("Terrain Slope [deg]")
    ax.set_zlabel("Maximum Velocity [m/s]")

    if save_plots:
        plt.savefig('plots/analysis_combined_terrain.png')

    plt.show()


if __name__ == '__main__':
    main()
