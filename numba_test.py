import subfunctions as sf
import nb_subfunctions as nb_sf
import define_experiment
import numpy as np
import define_rovers
import time
from scipy.interpolate import interp1d

def test_get_gear_ratio():
    # time the function

    speed_reducer = {
        'type': 'reverted',
        'diam_pinion': 0.1,
        'diam_gear': 0.2
    }

    nb_sf.get_gear_ratio(speed_reducer)

    tic = time.time()
    for i in range(100):
        sf.get_gear_ratio(speed_reducer)
    toc = time.time()
    sf_time = toc - tic

    tic = time.time()
    for i in range(100):
        nb_sf.get_gear_ratio(speed_reducer)
    toc = time.time()
    nb_time = toc - tic

    print(f'get_gear_ratio: sf_time = {sf_time}, nb_time = {nb_time}')

def test_tau_dc_motor():
    motor = define_rovers.rover1()['wheel_assembly']['motor']
    omega = np.linspace(0, 100, 1000, dtype=np.float64)

    nb_sf.tau_dcmotor(omega, motor)

    sf_time = None
    tic = time.time()
    for i in range(1000):
        sf.tau_dcmotor(omega, motor)
    toc = time.time()
    sf_time = toc - tic


    tic = time.time()
    for i in range(1000):
        nb_sf.tau_dcmotor(omega, motor)
    toc = time.time()
    nb_time = toc - tic

    print(f'tau_dcmotor: sf_time = {sf_time}, nb_time = {nb_time}')
    print(f'numba is {sf_time/nb_time} times faster')

def test_F_rolling():
    rover = define_rovers.rover1()
    planet = define_rovers.planet1()
    Crr = 0.3
    omega = np.linspace(0, 10, 100, dtype=np.float64)
    # pick random terrain angles -60 to 60 degrees the size of omega
    terrain_angle = np.random.rand(len(omega))*120 - 60

    nb_sf.F_rolling(omega, terrain_angle, rover, planet, Crr)

    tic = time.time()
    for i in range(1000):
        sf_answer = sf.F_rolling(omega, terrain_angle, rover, planet, Crr)
    toc = time.time()
    sf_time = toc - tic

    tic = time.time()
    for i in range(1000):
        nb_sf_answer = nb_sf.F_rolling(omega, terrain_angle, rover, planet, Crr)
    toc = time.time()
    nb_time = toc - tic

    # check that the answers are the same
    assert np.allclose(sf_answer, nb_sf_answer)

    print(f'F_rolling: sf_time = {sf_time}, nb_time = {nb_time}')
    print(f'numba is {sf_time/nb_time} times faster')

def test_rover_dynamics():
    rover = define_rovers.rover1()
    planet = define_rovers.planet1()
    experiment, endevent = define_experiment.experiment1()
    t = 20
    y = np.array([.25, 500])

    alpha_fun = interp1d(experiment['alpha_dist'], experiment['alpha_deg'], kind = 'cubic', fill_value = 'extrapolate') 

    tic = time.time()
    for i in range(1000):
        cache_answer = nb_sf.rover_dynamics(t, y, rover, planet, experiment, alpha_fun)
    toc = time.time()
    cache_time = toc - tic

    tic = time.time()
    for i in range(1000):
        no_cache_answer = nb_sf.rover_dynamics(t, y, rover, planet, experiment)
    toc = time.time()
    no_cache_time = toc - tic

    # check that the answers are the same
    assert np.allclose(cache_answer, no_cache_answer)
    print(f'rover_dynamics: cache_time = {cache_time}, no_cache_time = {no_cache_time}')
    print(f'cache is {no_cache_time/cache_time} times faster')

def test_simulate_rover():
    rover = define_rovers.rover1()
    planet = define_rovers.planet1()
    experiment, endevent = define_experiment.experiment1()

    quick_rover_dynamics = nb_sf.quick_dynamics_factory(rover, planet, experiment)

    tic = time.time()
    for i in range(10):
        nb_sf_answer = nb_sf.simulate_rover(rover, planet, experiment, endevent, quick_rover_dynamics=quick_rover_dynamics)
    toc = time.time()
    nb_time = toc - tic

    tic = time.time()
    for i in range(10):
        sf_answer = sf.simulate_rover(rover, planet, experiment, endevent)
    toc = time.time()
    sf_time = toc - tic

    # check that the answers are the same
    print(f'simulate_rover: sf_time = {sf_time}, nb_time = {nb_time}')
    print(f'optimization is {sf_time/nb_time} times faster')

def run_simulate_rover():
    rover = define_rovers.rover1()
    planet = define_rovers.planet1()
    experiment, endevent = define_experiment.experiment1()
    for i in range(50):
        nb_sf_answer = nb_sf.simulate_rover(rover, planet, experiment, endevent)

if __name__ == '__main__':
    # test_get_gear_ratio()
    # test_tau_dc_motor()
    # test_F_rolling()
    # test_rover_dynamics()
    run_simulate_rover()