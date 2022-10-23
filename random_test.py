from subfunctions import *
import define_rovers
import numpy as np

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

if __name__ == '__main__':
    # test_1()
    test_motorW()