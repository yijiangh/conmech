"""
Utility functions for reading/writing analysis result json files, and comparison.

"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import os
import json
import datetime
from collections import OrderedDict

import numpy as np

RES_UNIT_STANDARD = {
    "length_unit": "meter",
    "rot_angle_unit": "rad",
    "force_unit": "kN",
    "moment_unit": "kN-m",
}

def read_frame_analysis_result_json(res_file_path, verbose=False):
    """
    
    Parameters
    ----------
    res_file_path : [type]
        [description]
    verbose : bool, optional
        [description], by default False
    
    Returns
    -------
    [type]
        [description]
    """
    with open(res_file_path, 'r') as f:
        json_data = json.loads(f.read())
        if verbose: print('conmech result parsed: {}'.format(res_file_path))

    nD_data = json_data['node_displacement']
    nD = {}
    for nd in nD_data:
        # meter / rad
        nD[nd['node_id']] = nd['displacement']

    fR_data = json_data['fixity_reaction']
    fR = {}
    for fr in fR_data:
        fR[fr['node_id']] = np.array(fr['reaction'])

    eR = {}

    unit_info = {}
    for key_name in RES_UNIT_STANDARD.keys():
        unit_info[key_name] = json_data[key_name] if key_name in json_data else RES_UNIT_STANDARD[key_name]

    return nD, fR, eR, unit_info


def compare_analysis_results(res_A_path, res_B_path):
    pass