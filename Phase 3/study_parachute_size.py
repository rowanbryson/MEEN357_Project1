import numpy as np
import matplotlib.pyplot as plt
from define_edl_system import *
from subfunctions_EDL import *
from define_planet import *
from define_mission_events import *


def setup_edl(edl_system):
    # Overrides what might be in the loaded data to establish our desired
    # initial conditions
    edl_system['altitude'] = 11000    # [m] initial altitude
    edl_system['velocity'] = -578     # [m/s] initial velocity
    edl_system['rocket']['on'] = False
    edl_system['parachute']['deployed'] = True   # our parachute is open
    edl_system['parachute']['ejected'] = False   # and still attached
    edl_system['heat_shield']['ejected'] = False
    edl_system['sky_crane']['on'] = False
    edl_system['speed_control']['on'] = False
    edl_system['position_control']['on'] = False
    edl_system['rover']['on_ground'] = False # the rover has not yet landed
    
    # edl_system['parachute']['diameter'] = 17.0
    return edl_system

def get_specific_sim_runner(edl_system1, mars, mission_events, tmax, verbose):
    def sim_runner(parachute_diameter):
        edl_system2 = setup_edl(edl_system1)
        edl_system2['parachute']['diameter'] = parachute_diameter
        [T, Y, edl_system] = simulate_edl(edl_system2, mars, mission_events, tmax, verbose)

        end_time = T[-1]
        end_velocity = Y[0, -1]
        end_altitude = Y[1, -1]
        end_rover_rel_velocity = Y[5, -1]
        end_rover_rel_position = Y[6, -1]
        end_rover_position = end_altitude + end_rover_rel_position
        end_rover_velocity = end_velocity + end_rover_rel_velocity

        success = True
        # check if end_rover_position is within 1 m of the ground
        if abs(end_rover_position) > 1:
            success = False
        if end_altitude <= edl_system["sky_crane"]["danger_altitude"]:
            success = False
        if abs(end_rover_velocity) >= abs(edl_system["sky_crane"]["danger_speed"]):
            success = False
        return {'end_time': end_time, 'rover_end_speed': end_rover_velocity, 'success': success}

    return sim_runner

def main(show=True, save=False, verbose=False):
    # set up initial conditions
    edl_system = define_edl_system_1()
    mars = define_planet()
    mission_events = define_mission_events()
    tmax = 2000   # [s] maximum simulated time

    # set up a simulation runner that will run a simulation with a given parachute diameter
    sim_runner = get_specific_sim_runner(edl_system, mars, mission_events, tmax, False)

    parachute_diameters = list(np.arange(14, 19.5, 0.5))


    # run the simulations and record the results
    end_times, end_speeds, successes = [], [], []
    for parachute_diameter in parachute_diameters:
        end_info = sim_runner(parachute_diameter)
        end_times.append(end_info['end_time'])
        end_speeds.append(end_info['rover_end_speed'])
        successes.append(end_info['success'])

    # print the results
    if verbose:
        print('Parachute diameter [m], end time [s], end speed [m/s], success')
        for i in range(len(parachute_diameters)):
            print(f'{parachute_diameters[i]:.1f}, {end_times[i]:.1f}, {end_speeds[i]:.1f}, {successes[i]}')

    fig, ax = plt.subplots(3, 1)
    ax[0].plot(parachute_diameters, end_times)
    ax[0].set_title('Simulated time [s]')
    ax[0].set_ylabel('Time [s]')
    ax[1].plot(parachute_diameters, end_speeds)
    ax[1].set_title('Ending rover speed [m/s]')
    ax[1].set_ylabel('Speed [m/s]')

    # make a thin line to show the successes
    ax[2].plot(parachute_diameters, successes, color='black', linewidth=0.5)
    for parachute_diameter, success in zip(parachute_diameters, successes):
        if success:
            ax[2].plot(parachute_diameter, 1, 'go')
        else:
            ax[2].plot(parachute_diameter, 0, 'ro')

    ax[2].set_title('Rover landing success')
    ax[2].set_xlabel('Parachute diameter [m]')
    ax[2].set_ylabel('[1/0 = yes/no]')
    fig.suptitle('Parachute size study')
    plt.tight_layout()

    if save:
        plt.savefig('parachute_size_study.png')
    if show:
        plt.show()

if __name__ == '__main__':
    main(save=True)
