from subfunctions import *
import matplotlib.pyplot as plt
import numpy as np


def main():
    
    # Get Gear Ratio
    gear_ratio = get_gear_ratio(MARVIN_DICT['rover']['wheel_assembly']['speed_reducer'])
    print(gear_ratio)

    # Get Motor Torque
    motor_shaft_speed = np.linspace(0.0, 4.0, num=100)
    motor_torque = tau_dcmotor(motor_shaft_speed, MARVIN_DICT['rover']['wheel_assembly']['motor'])

    # Convert Motor Torque to Speed Reducer Torque
    sr_shaft_speed = motor_shaft_speed * gear_ratio
    sr_shaft_torque = motor_torque * gear_ratio
    sr_shaft_power = sr_shaft_speed * sr_shaft_torque
    
    #Plot 3 Graphs in 3x1 Grid
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1)

    ax1.plot(sr_shaft_torque, sr_shaft_speed)
    ax1.set_title('Speed Reducer Shaft Speed vs. Shaft Torque')
    ax1.set_ylabel('SR Shaft Speed (rad/s)')
    ax1.set_xlabel('SR Torque (N-m)')

    ax2.plot(sr_shaft_torque, sr_shaft_power)
    ax2.set_title('Speed Reducer Shaft Power vs. Shaft Torque')
    ax2.set_ylabel('SR Shaft Power (W)')
    ax2.set_xlabel('SR Torque (N-m)')

    ax3.plot(sr_shaft_speed, sr_shaft_power)
    ax3.set_title('Speed Reducer Shaft Power vs. Shaft Speed')
    ax3.set_ylabel('SR Shaft Power (W)')
    ax3.set_xlabel('SR Shaft Speed (rad/s)')

    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
        main()
