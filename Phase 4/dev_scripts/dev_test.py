import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from rich import print
import sampling_tools as st
import define_experiment as def_ex
import define_opt_problem as def_opt

def test_terrain_generator():
    terrain_designer = st.terrain_generator()
    terrain = next(terrain_designer)
    print(terrain)

def test_design_generator():
    CONTINUOUS_OPTIONS = def_opt.define_continuous_options()
    INTEGER_OPTIONS = def_opt.define_integer_options()
    CATEGORY_OPTIONS = def_opt.define_category_options()
    CONSTRAINTS = def_opt.define_constraints()
    designer = st.design_generator(CONTINUOUS_OPTIONS, INTEGER_OPTIONS, CATEGORY_OPTIONS, CONSTRAINTS)
    design = next(designer)
    print(design)

if __name__ == '__main__':
    # test_terrain_generator()
    # initialize the experiment
    test_design_generator()
