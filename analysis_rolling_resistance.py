from subfunctions import *
import matplotlib.pyplot as plt
import numpy as np


def main(save_plots=False):
    
    #Initialize Variables
    terrain_angle = 0
    omega = 2.0
    Crr_array = np.linspace(0.01, 0.4, num=25)
    forces = np.zeros(len(Crr_array))

    #Calculate Net Forces
    for i in range(len(Crr_array)):
        forces[i] = F_net(omega, terrain_angle, MARVIN_DICT['rover'], MARVIN_DICT['planet'], Crr_array[i])

    #Calculate Net Accelerations
    accelerations = forces / get_mass(MARVIN_DICT['rover'])

    #Plot acceleration vs. Crr
    plt.plot(Crr_array, accelerations)
    plt.title('Acceleration vs. Rolling Resistance Coefficient')
    plt.ylabel('Acceleration (m/s^2)')
    plt.xlabel('Crr')
    if save_plots:
        plt.savefig('plots/analysis_rolling_resistance.png')
    plt.show()


if __name__ == '__main__':
        main()