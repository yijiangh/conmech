import numpy as np
import conmech_py as cm
import os
cwd = os.getcwd()

json_path = cwd + "/sf-test_3-frame.json"

sc = cm.stiffness_checker(json_file_path = json_path, verbose = True)

ext_load = np.zeros([1,7])
ext_load[0,0] = 3
ext_load[0,3] = -500 * 1e3


sc.set_nodal_load(nodal_forces = ext_load, include_self_weight = False)
sc.solve()

# existing_ids = [0,1]
# sc.solve(existing_ids)