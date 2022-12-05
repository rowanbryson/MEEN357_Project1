import sampling_tools as samp_tools
import simulation_tools as sim_tools
import define_experiment as def_ex
import define_opt_problem as def_opt
import json
from rich import print

import numpy as np
import scipy
import scipy.optimize as opt

from loguru import logger

LOG_FILEPATH = 'Phase 4/opt_logs.txt'

logger.add(sink=LOG_FILEPATH, format='{time:YYYY-MM-DD HH:mm:ss} | {level} | {message}', level='INFO')

# initialize the experiment
CONTINUOUS_OPTIONS = def_opt.define_continuous_options()
INTEGER_OPTIONS = def_opt.define_integer_options()
CATEGORY_OPTIONS = def_opt.define_category_options()
CONSTRAINTS = def_opt.define_constraints()

def init_empty_design_file():
    with open('Phase 4/top_10_designs.json', 'w') as f:
        json.dump([], f)


def submit_design_candidate(nums, design, time, filepath):
    """
    opens the json file at filepath and submits the design candidate

    if it is one of the top 10 designs, it adds to the list
    if there were already 10 designs, it also removes the worst one
    """
    with open(filepath, 'r') as f:
        # check if the design qualifies as top 10
        data = json.load(f)
        if len(data) < 10:
            data.append({
                'design_nums': list(nums),
                'design': design,
                'time': time,
            })
        else:
            # check if the design is better than any of the current designs
            for i, design_data in enumerate(data):
                if time < design_data['time']:
                    data[i] = {
                        'design_nums': nums,
                        'design': design,
                        'time': time,
                    }
                    break
            else:
                # if the design is not better than any of the current designs, don't add it
                return
            # sort the designs by time
            data = sorted(data, key=lambda x: x['time'])

    # write the new data to the file
    with open(filepath, 'w') as f:
        json.dump(data, f, indent=4)


def run_optimizer(popsize, maxiter):
    experiment_holder = samp_tools.ExperimentHolder()

    # use scipy differential evolution
    bounds = opt.Bounds([0] * 9, [1] * 9)
    # evaluator = sim_tools.simplified_evaluator_factory(
    #     experiment, end_event, CONTINUOUS_OPTIONS, INTEGER_OPTIONS, CATEGORY_OPTIONS, CONSTRAINTS, 1000, False)
    evaluator = sim_tools.multi_terrain_evaluator_factory(experiment_holder, CONTINUOUS_OPTIONS, INTEGER_OPTIONS, CATEGORY_OPTIONS, CONSTRAINTS, 1000, False) 

    design_generator = samp_tools.design_generator(
        CONTINUOUS_OPTIONS, INTEGER_OPTIONS, CATEGORY_OPTIONS, CONSTRAINTS)
    x0, design = next(design_generator)

    # set up initial population
    designer = samp_tools.design_generator(
        CONTINUOUS_OPTIONS, INTEGER_OPTIONS, CATEGORY_OPTIONS, CONSTRAINTS)
    initial_design_nums = np.vstack(
        [next(designer)[0] for i in range(popsize)])

    def callback(xk, convergence):
        experiment_holder.update_experiments(5)
        design = samp_tools.nums_to_design(
            xk, CONTINUOUS_OPTIONS, INTEGER_OPTIONS, CATEGORY_OPTIONS)
        print('design:', design)
        print('convergence:', convergence)
        print()

    print('[green]Starting Optimizer...[/]')
    result = opt.differential_evolution(evaluator, bounds, maxiter=maxiter, popsize=popsize,
                                        tol=0.01, polish=False, disp=True, init=initial_design_nums, callback=callback)

    design_nums = result.x
    design = samp_tools.nums_to_design(
        design_nums, CONTINUOUS_OPTIONS, INTEGER_OPTIONS, CATEGORY_OPTIONS)

    message = result.message
    success = result.success

    

    # evaluate the design on many terrains to guage its performance
    experiment_holder.update_experiments(20)
    successes = 0
    total_time = 0
    for experiment, end_event in experiment_holder.experiments:
        result = sim_tools.evaluate_design(design, experiment, end_event, CONTINUOUS_OPTIONS, INTEGER_OPTIONS, CATEGORY_OPTIONS, CONSTRAINTS, 1000, full_output=True)
        if result['success']:
            successes += 1
        total_time += result['time']
    percent_success = successes / 20
    avg_time = total_time / 20
        

    print(f'scipy optimizer message: {message}')
    print(f'optimizer success: {success}')
    print('design', design)
    print('design_nums', design_nums)

    print(f'percent success: {percent_success}')
    print(f'avg time: {avg_time}')

    info_to_log = {
        'time': avg_time,
        'percent_success': percent_success,
        'design_nums': design_nums,
        'design': design,
        'optimizer_message': message,
    }

    logger.info(info_to_log)
    submit_design_candidate(design_nums, design,
                            avg_time, 'Phase 4/top_10_designs.json')


if __name__ == '__main__':
    logger.info('Starting optimization session...')
    for i in range(10):
        try:
            run_optimizer(50, 50)
        except Exception as e:
            logger.error(e)