# allow imports from parent directory
import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import subfunctions as sf
import numpy as np
import matplotlib.pyplot as plt
import define_rovers, define_experiment
from rich import print
import experiment_visualization

def get_telemetry_data():
    '''initialize the rover and experiment, then run the experiment with simulate_rover'''
    rover = define_rovers.rover1()
    experiment, default_end_event = define_experiment.experiment1()
    planet = define_rovers.planet1()
    end_event = {
        'max_distance' : 1000,
        'max_time' : 10000,
        'min_velocity' : 0.01
    }
    rover = sf.simulate_rover(rover, planet, experiment, end_event)
    return rover['telemetry']

def plot_data(telemetry, show_plots=True, save_plots=False, plot_terrain=False):
    '''plot the telemetry data and show the plot
    
    Parameters
    ----------
    telemetry : dict
        telemetry data from the rover
    show_plots : bool, optional
        whether to show the plots, by default True
    save_plots : bool, optional
        whether to save the plots, by default False
    plot_terrain : bool, optional
        whether to plot the terrain, by default False
        ! this is useful for debugging, but not required for grading, so it is False by default
    '''
    time = telemetry['time']
    fig1, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True)

    ax1.plot(time, telemetry['position'])
    ax1.set_ylabel('Position [m]')
    ax2.plot(time, telemetry['velocity'])
    ax2.set_ylabel('Velocity [m/s]')
    ax3.plot(time, telemetry['power'])
    ax3.set_ylabel('Power [W]')
    ax3.set_xlabel('Time [s]')
    fig1.suptitle('Rover Telemetry Graphs')
    ax1.grid()
    ax2.grid()
    ax3.grid()

    # also plot the terrain information if specified
    if plot_terrain:
        ax = experiment_visualization.height_over_x_distance(define_experiment.experiment1()[0])
        plt.sca(ax)
        fig2 = plt.gcf()
    if show_plots:
        plt.show()
    if save_plots:
        fig1.savefig('plots/phase_2/telemetry_graphs.png')
        if plot_terrain:
            fig2.savefig('plots/phase_2/terrain_graph.png')


def save_telemetry_data(telemetry, filepath='plots/phase_2/telemetry_data.csv'):
    '''save the telemetry data to a csv file

    Parameters
    ----------
    telemetry : dict
        telemetry data from the rover
    filepath : str, optional
        path to save the csv file, by default 'plots/phase_2/telemetry_data.csv'
    '''
    fields_to_save = ['completion_time', 'distance_traveled', 'max_velocity', 'average_velocity', 'battery_energy', 'energy_per_distance']
    # save a csv file with the telemetry data
    with open(filepath, 'w') as f:
        for field in fields_to_save:
            f.write(f'{field},{telemetry[field]}\n')

if __name__ == '__main__':
    telemetry = get_telemetry_data()
    plot_data(telemetry, show_plots=True, save_plots=True, plot_terrain=True)
    # save_telemetry_data(telemetry)
