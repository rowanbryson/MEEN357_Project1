from subfunctions import *
from define_rovers import rover1
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

rover = rover1()

# Assign variables to effcy and effcy_tau
effcy_tau_values = rover['wheel_assembly']['motor']['effcy_tau']
effcy_values = rover['wheel_assembly']['motor']['effcy']

# Interpolate the efficiency data
effcy_fun = interp1d(effcy_tau_values, effcy_values, kind='cubic')

# Plot the interpolated efficiency data (100 points, star shaped markers)
tau = np.linspace(0, 165, 100)
plt.plot(tau, effcy_fun(tau), 'r*-')
plt.title('Interpolated Efficiency Data')
plt.xlabel('Torque [Nm]')
plt.ylabel('Efficiency [-]')
plt.grid()
plt.savefig('efficiency_visualization.png')
plt.show()