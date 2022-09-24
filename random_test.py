from subfunctions import *
import numpy as np



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