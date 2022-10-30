from turtle import title
import define_experiment
from subfunctions import *
import define_rovers
import numpy as np
from functools import partial
import warnings
try:
    from rich import print
except:
    warnings.warn('I am using the print function from the rich module to make the output look nicer. You can install it with "pip install rich" -Jae')


# this file is for testing functions quickly without writing test cases

def test_1():
    MARVIN_DICT = {'rover': define_rovers.rover1(), 'planet': define_rovers.planet1()}
    omega = np.array([0.0, 0.5, 1.0, 2.0, 3.0, 3.8])
    terrain = np.array([-5.0, 0.0, 5.0, 10.0, 20.0, 30.0])
    rover = MARVIN_DICT['rover']
    planet = MARVIN_DICT['planet']
    Crr = 0.1
    print(F_drive(omega, rover))
    print(F_gravity(terrain, rover, planet))
    print(F_rolling(omega, terrain, rover, planet, Crr))
    print(F_net(omega, terrain, rover, planet, Crr))

    # Output:
    # https://imgur.com/a/pfgCgkn

def test_motorW():
    rover = define_rovers.rover1()

    gear_ratio = get_gear_ratio(rover['wheel_assembly']['speed_reducer'])
    wheel_radius = rover['wheel_assembly']['wheel']['radius']

    v = np.array([0, 1, 2, 3, 4, 5, 6, 7])

    print(f'gear_ratio: {gear_ratio}')
    print(f'wheel_radius: {wheel_radius}')
    print(f'v: {v}')
    print(f'motorW(v): {motorW(v, rover)}')

def test_simulate_rover():
    rover = define_rovers.rover1()
    planet = define_rovers.planet1()
    experiment, end_event = define_experiment.experiment1()
    rover_after = simulate_rover(rover, planet, experiment, end_event)
    print(rover_after['telemetry'])

def test_rover_dynamics():
    rover = define_rovers.rover1()
    planet = define_rovers.planet1()
    experiment, end_event = define_experiment.experiment2()
    rover_dp = partial(rover_dynamics, rover=rover, planet=planet, experiment=experiment)

    print(rover_dp(0, np.array([0, 0])))

def test_batt_energy():
    rover = define_rovers.rover1()
    t = np.array([0, 1, 2, 3, 4, 5, 6])
    v = np.array([0.33, 0.32, 0.33, 0.2, 0.2, 0.25, 0.28])
    print(battenergy(t, v, rover, quick_plot=True))

def test_mech_power():
    rover = define_rovers.rover1()
    v = np.linspace(0, 10, 100)
    plt.plot(v)
    plt.ylabel('velocity [m/s]')
    print(mechpower(v, rover, quick_plot=True))

def test_taudcmotor():
    rover = define_rovers.rover1()
    motor = rover['wheel_assembly']['motor']
    omega = np.linspace(-1, 4, 100)
    torque = tau_dcmotor(omega, motor)
    plt.plot(omega, torque)
    plt.xlabel('omega [rad/s]')
    plt.ylabel('torque [Nm]')
    plt.show()

def get_motor_effcy():
    v = 0.25
    rover = define_rovers.rover1()
    motor = rover['wheel_assembly']['motor']
    motor_effcy_interp = interp1d(motor['effcy_tau'], motor['effcy'], kind='cubic')
    print(motor_effcy_interp(v))

if __name__ == '__main__':
    # test_1()
    # test_motorW()
    test_simulate_rover()
    # test_rover_dynamics()
    # test_batt_energy()
    # test_mech_power()
    # test_taudcmotor()
    # get_motor_effcy()