#!/usr/bin/env python
from __future__ import print_function

import os
import argparse
import numpy as np
from numpy.linalg import norm
from termcolor import cprint
import random
import time

import matplotlib
import matplotlib.pyplot as plt
from collections import OrderedDict, defaultdict

from pyconmech import StiffnessChecker

# TEST_JSON_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'tests', 'test_data')
TEST_JSON_PATH = 'C:/Users/yijiangh/Documents/pb_ws/pb-construction/assembly_instances/extrusion'
EXCLUDE_LIST = []

##################################################

#ALPHA = None
#EDGES = ['face', 'face', 'g', 'b']
#COLORS = ['r', 'y', 'none', 'none']

MARKERS = ['x', 'o', '+']
EDGES = ['face', 'g', 'face', ] # C2
#COLORS = ['C0', 'C1', 'none']
#COLORS = ['c', 'y', 'none'] # m
COLORS = ['r', 'none', 'b']

# https://matplotlib.org/api/markers_api.html
# https://matplotlib.org/2.0.2/api/colors_api.html

def scatter_plot(data):
    all_sizes = sorted({result['nE'] for _, result in data.items()})
    print('Sizes:', all_sizes)

    sizes = []
    np_runtimes = []
    cpp_runtimes = []
    np_init_times = []
    cpp_init_times = []
    for _, d in data.items():
        sizes.append(d['nE'])
        np_init_times.append(d['np'][0])
        np_runtimes.append(d['np'][1])
        cpp_init_times.append(d['cpp'][0])
        cpp_runtimes.append(d['cpp'][1])

    max_time = np.max(np.hstack([np_init_times, cpp_init_times, np_runtimes, cpp_runtimes]))

    fig, axes = plt.subplots(1, 2)

    # runtime plot
    # ticks
    buffer_time = 0.01*max_time

    axes[0].scatter(all_sizes, np.zeros(len(all_sizes)), marker='|', color='k') # black
    axes[0].scatter(sizes, np_runtimes, marker=MARKERS[0],
                color=COLORS[0], edgecolors=EDGES[0],
                alpha=0.75, label='numpy')
    axes[0].scatter(sizes, cpp_runtimes, marker=MARKERS[1],
                color=COLORS[1], edgecolors=EDGES[1],
                alpha=0.75, label='cpp')
    axes[0].set_title('solve time')
    #axes[0].xticks(range(1, max_size+1)) #, [get_index(problem) for problem in problems])
    axes[0].set_xlim([1, 1000]) #max(all_sizes)])
    # axes[0].set_ylim([0, max_time+buffer_time]) #max(all_sizes)])
    axes[0].set_xlabel('# elements')
    axes[0].set_ylabel('runtime (sec)')
    axes[0].legend()
    #plt.legend(loc='upper left')
    # plt.legend(loc='upper center')
    #plt.savefig('test')

    # init time plot
    axes[1].scatter(all_sizes, np.zeros(len(all_sizes)), marker='|', color='k') # black
    axes[1].scatter(sizes, np_init_times, marker=MARKERS[0],
                color=COLORS[0], edgecolors=EDGES[0],
                alpha=0.75, label='numpy')
    axes[1].scatter(sizes, cpp_init_times, marker=MARKERS[1],
                color=COLORS[1], edgecolors=EDGES[1],
                alpha=0.75, label='cpp')
    axes[1].set_title('Initialization time')
    axes[1].set_xlim([1, 1000]) #max(all_sizes)])
    # axes[1].set_ylim([0, max_time+buffer_time]) #max(all_sizes)])
    axes[1].set_xlabel('# elements')
    axes[1].set_ylabel('runtime (sec)')
    axes[1].legend()

    plt.show()
    # logarithmic scale

##################################################

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

    nE = len(np_sc.elements)
    full_eids = list(range(nE))

    for i in range(n_attempts):
        # nE_exist = random.sample(full_eids, 1)[0]
        # existing_ids = random.sample(full_eids, nE_exist)
        existing_ids = []
        if verbose: 
            print('-----\nIter {} : exist E#{}/{}'.format(i, len(existing_ids), nE))

        cpp_start_time = time.time()
        cpp_success = cpp_sc.solve(existing_ids, if_cond_num=True)
        cpp_t_tol, _ = cpp_sc.get_nodal_deformation_tol()
        if cpp_success:
            cpp_max_t, _, cpp_max_tid, _ = cpp_sc.get_max_nodal_deformation()
        accum_cpp_time += (time.time() - cpp_start_time)

        np_start_time = time.time()
        np_success = np_sc.solve(existing_ids, if_cond_num=True)
        np_t_tol, _ = np_sc.get_nodal_deformation_tol()
        if np_success:
            np_max_t, _, np_max_tid, _ = np_sc.get_max_nodal_deformation()
        accum_np_time += (time.time() - np_start_time)

        if verbose: 
            if cpp_success:
                cprint('{} | Success: {} | Max trans deformation: Node #{}: {}|{}'.format('cpp', cpp_success, cpp_max_tid, cpp_max_t, cpp_t_tol), 'magenta')
            if np_success:
                cprint('{} | Success: {} | Max trans deformation: Node #{}: {}|{}'.format('numpy', np_success, np_max_tid, np_max_t, np_t_tol), 'cyan')
        assert np_success == cpp_success, 'np {} | cpp {}'.format(np_success, cpp_success)
    
    print('#'*15)
    print('File: {} | nE : {}'.format(np_sc.model_name, nE))
    cprint('Init time: \nnumpy: {} | cpp : {}'.format(np_init_time, cpp_init_time))
    cprint('Avg solve time (#{} run): \nnumpy: {} | cpp : {}'.format(n_attempts, accum_np_time/n_attempts, accum_cpp_time/n_attempts))
    data = {}
    data['model'] = np_sc.model_name
    data['nE'] = nE
    data['np'] = (np_init_time, accum_np_time/n_attempts)
    data['cpp'] = (cpp_init_time, accum_cpp_time/n_attempts)
    return data

def main():
    np.set_printoptions(precision=3)
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--n_trails', default=10, help='number of trails')
    parser.add_argument('-d', '--debug', action='store_true', help='Debug verbose mode')
    parser.add_argument('-p', '--problem', default='topopt-100.json', help='The json file name of the problem to solve')
    parser.add_argument('-v', '--verbose', action='store_true', help='verbose output')
    parser.add_argument('-a', '--all', action='store_true', help='check all instances in a folder')
    parser.add_argument('-w', '--write', action='store_true', help='generate statistic graph')
    args = parser.parse_args()
    print('Arguments:', args)

    stat_data = {}
    if not args.all:
        frame_file_path = os.path.join(TEST_JSON_PATH, args.problem)
        stat_data[args.problem] = compare_stiffness_engine(frame_file_path, n_attempts=int(args.n_trails), verbose=args.verbose)
    else:
        for filename in sorted(os.listdir(TEST_JSON_PATH)):
            if filename.endswith('.json') and 'load_case' not in filename and filename not in EXCLUDE_LIST:
                frame_file_path = os.path.join(TEST_JSON_PATH, filename)
                stat_data[filename] = compare_stiffness_engine(frame_file_path, n_attempts=int(args.n_trails), verbose=args.verbose)

    if args.write:
        # generate stat plot
        scatter_plot(stat_data)


if __name__ == '__main__':
    main()