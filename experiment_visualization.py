from subfunctions import *
from define_experiment import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.integrate import solve_ivp

def angle_over_distance(experiment):
    '''Plots the angle of the terrain over time
    
    This is the on we will be graded on.
    '''
    #initalization of variables
    alpha_dist = experiment['alpha_dist']
    alpha_deg = experiment['alpha_deg']

    #define the alpha function to find terrain angles
    alpha_fun = interp1d(alpha_dist, alpha_deg, kind='cubic', fill_value='extrapolate')

    #plot of terrain angle vs rover position
    xaxis = np.linspace(alpha_dist[0], alpha_dist[-1], 100)
    fig, ax = plt.subplots()
    ax.plot(xaxis, alpha_fun(xaxis), label="100 evaluated terrain angles")
    ax.plot(alpha_dist, alpha_deg, 'r*', label="define_experiment data")
    ax.legend()
    ax.set_title("Rover's Terrain Angle [deg] vs Position [m]")
    ax.set_xlabel("Position [m]")
    ax.set_ylabel("Terrain Angle (deg)")
    return ax

def height_over_x_distance(experiment):
    alpha_fun = interp1d(experiment['alpha_dist'], experiment['alpha_deg'], kind='cubic', fill_value='extrapolate')
    def terrain_derivatives(s, z):
        angle = alpha_fun(s)  # angle of terrain in degrees
        dxds = np.cos(np.deg2rad(angle))  # derivative of x wrt s
        dyds = np.sin(np.deg2rad(angle))  # derivative of y wrt s
        return np.array([dxds, dyds])
    
    t_span = [experiment['alpha_dist'][0], experiment['alpha_dist'][-1]]
    origin = np.array([0, 0])
    sol = solve_ivp(terrain_derivatives, t_span, y0=origin, method='RK45', t_eval=np.linspace(t_span[0], t_span[1], 10000))

    fig, ax = plt.subplots()
    ax.plot(sol.y[0], sol.y[1])
    ax.set_title("Terrain Height vs x Distance")
    ax.set_xlabel("x Distance [m]")
    ax.set_ylabel("Terrain Height [m]")
    return ax


def visualize(experiment, full=False):
    ax = angle_over_distance(experiment)
    plt.sca(ax)
    if full:
        ax = height_over_x_distance(experiment)
        plt.sca(ax)
        plt.show()

if __name__ == '__main__':
    experiment, end_event = experiment2()
    visualize(experiment, full=True)