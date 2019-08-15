from __future__ import print_function
import os

import pytest
import numpy as np
import random
from numpy.testing import assert_equal, assert_almost_equal
from pyconmech import stiffness_checker

def test_stiffness_checker_single_full_solve():
    cwd = os.getcwd()
    root_dir = os.path.dirname(os.path.abspath(__file__))
    json_path = os.path.join(root_dir, 'test_data', 'tower_3D.json')
    sc = stiffness_checker(json_file_path=json_path)

    sc.set_self_weight_load(True)
    assert sc.solve()
    # existing_ids = [0, 24, 25, 26, 27]
    # sc.solve(existing_ids)
