#!/usr/bin/env python
from __future__ import print_function

import os
import argparse
import numpy as np
from numpy.linalg import norm
from termcolor import cprint
import random
import time

from pyconmech import StiffnessChecker

TEST_JSON_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'tests', 'test_data')

def compare_stiffness_engine(frame_file_path, n_attempts=10, verbose=True):
    trans_tol = 1e-3
    parsed_pt_loads= None
    include_sw = True

    accum_np_time = 0
    accum_cpp_time = 0

    np_init_st_time = time.time()
    np_sc = StiffnessChecker.from_json(json_file_path=frame_file_path, checker_engine='numpy', verbose=False)
    np_sc.set_nodal_displacement_tol(trans_tol=trans_tol)
    np_init_time = (time.time() - np_init_st_time)
    # np_sc.set_loads(point_loads=parsed_pt_loads, include_self_weight=include_sw)

    cpp_init_st_time = time.time()
    cpp_sc = StiffnessChecker.from_json(json_file_path=frame_file_path, checker_engine='cpp', verbose=False)
    cpp_sc.set_nodal_displacement_tol(trans_tol=trans_tol)
    # cpp_sc.set_loads(point_loads=parsed_pt_loads, include_self_weight=include_sw)
    cpp_init_time = (time.time() - cpp_init_st_time)

    full_eids = list(range(len(np_sc.elements)))

    for _ in range(n_attempts):
        print('-----')
        nE_exist = random.sample(full_eids, 1)[0]
        # existing_ids = random.sample(full_eids, nE_exist)
        existing_ids = []
        print(existing_ids)
        
        np_start_time = time.time()
        np_success = np_sc.solve(existing_ids, if_cond_num=True)
        np_t_tol, _ = np_sc.get_nodal_deformation_tol()
        np_max_t, _, np_max_tid, _ = np_sc.get_max_nodal_deformation()
        accum_np_time += (time.time() - np_start_time)

        cprint('{} | Max trans deformation: Node #{}: {}|{}'.format('numpy', np_max_tid, np_max_t, np_t_tol), 'cyan')

        cpp_start_time = time.time()
        cpp_success = cpp_sc.solve(existing_ids, if_cond_num=True)
        cpp_t_tol, _ = cpp_sc.get_nodal_deformation_tol()
        cpp_max_t, _, cpp_max_tid, _ = cpp_sc.get_max_nodal_deformation()
        accum_cpp_time += (time.time() - cpp_start_time)

        cprint('{} | Max trans deformation: Node #{}: {}|{}'.format('cpp', cpp_max_tid, cpp_max_t, cpp_t_tol), 'magenta')

        assert np_success == cpp_success
    
    print('#'*15)
    cprint('Init time: \nnumpy: {} | cpp : {}'.format(np_init_time, cpp_init_time))
    cprint('Avg solve time: \nnumpy: {} | cpp : {}'.format(accum_np_time/n_attempts, accum_cpp_time/n_attempts))

def main():
    np.set_printoptions(precision=3)
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--n_trails', default=10, help='number of trails')
    parser.add_argument('-d', '--debug', action='store_true', help='Debug verbose mode')
    parser.add_argument('-p', '--problem', default='topopt-100.json', help='The json file name of the problem to solve')
    # parser.add_argument('-v', '--viewer', action='store_true', help='')
    args = parser.parse_args()
    print('Arguments:', args)

    frame_file_path = os.path.join(TEST_JSON_PATH, args.problem)
    compare_stiffness_engine(frame_file_path, n_attempts=int(args.n_trails), verbose=True)


if __name__ == '__main__':
    main()