import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
import define_experiment as def_ex
import subfunctions_Phase4 as sf

sys.path.append(os.path.join(os.path.dirname(__file__), '../..'))
from graphers import experiment_visualization as ev


experiment, _ = def_ex.experiment1()
print('experiment:', experiment)

ev.visualize(experiment, full=True)