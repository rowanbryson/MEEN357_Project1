import subfunctions as sf
import define_experiment
import define_rovers


if __name__ == '__main__':
    rover = define_rovers.rover1()
    planet = define_rovers.planet1()
    experiment, endevent = define_experiment.experiment1()
    for i in range(50):
        sf.simulate_rover(rover, planet, experiment, endevent)