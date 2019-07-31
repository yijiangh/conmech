# should be refactored...
# better way to include the input json file for tests?
# pytest

import numpy as np
import pyconmech as cm
import os
cwd = os.getcwd()

file_dir =\
"/Users/yijiangh/Dropbox (MIT)/Projects/conmech/conmech/src/stiffness_checker/matlab/test/problem_instances/";

file_name = "tower_3D_wpillar.json";
file_path = file_dir + file_name;

sc = cm.stiffness_checker(json_file_path = file_path, verbose = True)

sc.set_self_weight_load(True)
sc.set_nodal_displacement_tol(transl_tol = 1e-3, rot_tol = 3*(3.14/360))

sc.solve()

# existing_ids = [0, 24, 25, 26, 27]
# sc.solve(existing_ids)
