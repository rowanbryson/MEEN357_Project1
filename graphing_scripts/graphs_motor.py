import sys
sys.path.append('../')
import define_rovers
from subfunctions import *
import matplotlib.pyplot as plt
import numpy as np

# this is a band aid to make the code work, the new subfunctions.py file doesn't define MARVIN_DICT
# we might eventually want to change this to a function that takes in a rover name
MARVIN_DICT = define_rovers.default_data_dict

def main(save_plots=False):
    motor_shaft_speed = np.linspace(0.0, 4.0, num=100)
    motor_shaft_torque = tau_dcmotor(motor_shaft_speed, MARVIN_DICT['rover']['wheel_assembly']['motor'])
    motor_power = motor_shaft_speed*motor_shaft_torque

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1)

    ax1.plot(motor_shaft_torque, motor_shaft_speed)
    ax1.set_title("Motor shaft torque [Nm] vs Motor shaft speed [rad/s]")
    ax1.set_xlabel("Motor shaft torque [Nm]")
    ax1.set_ylabel("Motor shaft speed [rad/s]")

    ax2.plot(motor_shaft_torque, motor_power)
    ax2.set_title("Motor shaft torque [Nm] vs Motor power [W]")
    ax2.set_xlabel("Motor shaft torque [Nm]")
    ax2.set_ylabel("Motor power [W]")

    ax3.plot(motor_shaft_speed, motor_power)
    ax3.set_title("Motor shaft speed [rad/s] vs Motor power [W]")
    ax3.set_xlabel("Motor shaft speed [rad/s]")
    ax3.set_ylabel("Motor power [W]")

    plt.tight_layout()
    if save_plots:
        plt.savefig('plots/graphs_motor.png')
    plt.show()

if __name__ == '__main__':
        main()