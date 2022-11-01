import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from subfunctions import *
from define_experiment import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.integrate import solve_ivp

def angle_over_distance(experiment, ax=None):
    '''Plots the angle of the terrain over time
    
    This is the on we will be graded on.

    Inputs:
    -------
    experiment : dict
        The experiment dictionary
    ax : matplotlib.axes.Axes, optional
        The axes to plot on. If None, a new figure is created.
    '''
    #initalization of variables
    alpha_dist = experiment['alpha_dist']
    alpha_deg = experiment['alpha_deg']

    #define the alpha function to find terrain angles
    alpha_fun = interp1d(alpha_dist, alpha_deg, kind='cubic', fill_value='extrapolate')

    #plot of terrain angle vs rover position
    xaxis = np.linspace(alpha_dist[0], alpha_dist[-1], 100)
    if ax is None:
        fig, ax = plt.subplots()
    ax.plot(xaxis, alpha_fun(xaxis), label="100 evaluated terrain angles")
    ax.plot(alpha_dist, alpha_deg, 'r*', label="define_experiment data")
    ax.legend()
    ax.set_title("Rover's Terrain Angle [deg] vs Position [m]")
    ax.set_xlabel("Position [m]")
    ax.set_ylabel("Terrain Angle (deg)")
    return ax

def height_over_x_distance(experiment, ax=None):
    '''
    Plots the height of the terrain with respect to the rover's position on the x axis

    This one is not required, but it is a nice visualization of the terrain,
    and it is helpful for deciding wether our output makes sense.
    '''
    alpha_fun = interp1d(experiment['alpha_dist'], experiment['alpha_deg'], kind='cubic', fill_value='extrapolate')
    def terrain_derivatives(s, z):
        angle = alpha_fun(s)  # angle of terrain in degrees
        dxds = np.cos(np.deg2rad(angle))  # derivative of x wrt s
        dyds = np.sin(np.deg2rad(angle))  # derivative of y wrt s
        return np.array([dxds, dyds])
    
    t_span = [experiment['alpha_dist'][0], experiment['alpha_dist'][-1]]
    origin = np.array([0, 0])
    sol = solve_ivp(terrain_derivatives, t_span, y0=origin, method='RK45', t_eval=np.linspace(t_span[0], t_span[1], 10000))

    if ax is None:
        fig, ax = plt.subplots()
    ax.plot(sol.y[0], sol.y[1])
    ax.set_title("Terrain Height vs x Distance")
    ax.set_xlabel("x Distance [m]")
    ax.set_ylabel("Terrain Height [m]")
    ax.grid()
    ax.set_aspect('equal', adjustable='box')
    return ax


def visualize(experiment, full=False, save_plot=False):
    '''
    This function handles calling the plot functions. It allows the
    plotting to be modular, which I needed because I wanted to be able to access
    the terrain ploting function from the rover simulation in order to plot the
    terrain right next to the telemetry data.

    Inputs:
    -------
    experiment : dict
        The experiment dictionary
    full : bool, optional
        If True, all plots are shown. If False, only the required plot (terrain angle) is shown.
    save_plot : bool, optional
        If True, the plot is saved to a file.
    '''
    ax = angle_over_distance(experiment)
    plt.sca(ax)
    if save_plot:
        plt.savefig('plots/phase_2/angle_over_distance.png')
    if full:
        ax = height_over_x_distance(experiment)
        plt.sca(ax)
        if save_plot:
            plt.savefig('plots/phase_2/height_over_x_distance.png')
    plt.show()

if __name__ == '__main__':
    experiment, end_event = experiment1()
    visualize(experiment, full=True, save_plot=True)