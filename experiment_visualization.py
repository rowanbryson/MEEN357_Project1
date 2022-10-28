from subfunctions import *
from define_experiment import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

experiment, end_event = experiment2()

#initalization of variables
alpha_dist = experiment['alpha_dist']
alpha_deg = experiment['alpha_deg']

#define the alpha function to find terrain angles
alpha_fun = interp1d(alpha_dist, alpha_deg, kind='cubic', fill_value='extrapolate')

#plot of terrain angle vs rover position
xaxis = np.linspace(alpha_dist[0], alpha_dist[-1], 100)
plt.plot(xaxis, alpha_fun(xaxis), label="100 evaluated terrain angles")
plt.plot(alpha_dist, alpha_deg, 'r*', label="define_experiment data")
plt.legend()
plt.title("Rover's Terrain Angle [deg] vs Position [m]")
plt.xlabel("Position [m]")
plt.ylabel("Terrain Angle (deg)")
plt.show()
