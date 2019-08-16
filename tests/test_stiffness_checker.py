from __future__ import print_function
import os
from tempfile import TemporaryDirectory

import pytest
import numpy as np
import random
from numpy.testing import assert_equal, assert_almost_equal

from pyconmech import stiffness_checker
from pyconmech import read_frame_json, write_frame_json
from pyconmech.database import MATERIAL_PROPERTY_DATABASE


def test_get_material_from_database():
    PLA_MAT = MATERIAL_PROPERTY_DATABASE['PLA']
    STEEL_S235_MAT = MATERIAL_PROPERTY_DATABASE['Steel-S235']
    print(PLA_MAT)
    print(STEEL_S235_MAT)


@pytest.mark.wip
def test_stiffness_checker_single_full_solve():
    file_name = 'tower_3D.json'
    root_dir = os.path.dirname(os.path.abspath(__file__))
    json_path = os.path.join(root_dir, 'test_data', file_name)
    sc = stiffness_checker(json_file_path=json_path)
    sc.set_self_weight_load(True)

    existing_ids = [0, 4, 7, 8, 9] # some parts are floating
    assert not sc.solve(existing_ids, eid_sanity_check=True)
    success, fail_nD, fail_fR, fail_eR = sc.get_solved_results() # can get, but result not useful

    assert sc.solve()
    success, nD, fR, eR = sc.get_solved_results()
    assert len(nD) == len(sc.node_points)
    assert len(eR) == len(sc.elements)
    assert len(fR) == len(sc.fix_node_ids)

    for fnid in sc.fix_node_ids:
        assert_equal(nD[fnid], np.zeros(6))


def test_init_stiffness_checker():
    file_name = 'tower_3D_broken_lines.json'

    here = os.path.dirname(os.path.abspath(__file__))
    file_path = os.path.join(here, 'test_data', file_name)

    print('init stiffness_checker from json file path...')
    # this is the same as: sc_from_json = stiffness_checker(file_path)
    sc_from_json = stiffness_checker.from_json(file_path, verbose=False)

    print('init stiffness_checker from frame data...')
    nodes, elements, fixed_node_ids, fix_specs, model_type, material_dict, model_name = \
        read_frame_json(file_path, verbose=True)
    sc_from_data = stiffness_checker.from_frame_data(nodes, elements, fixed_node_ids, material_dict, 
        fixity_specs=fix_specs, unit='meter', model_type=model_type, model_name=model_name, verbose=False)

    assert sc_from_json.model_name == sc_from_data.model_name
    assert sc_from_json.model_type == sc_from_data.model_type
    assert sc_from_json.material == sc_from_json.material
    for n1, n2 in zip(sc_from_json.node_points, sc_from_data.node_points):
        assert_equal(n1, n2)
    for e1, e2 in zip(sc_from_json.elements, sc_from_data.elements):
        assert_equal(e1, e2)
    for fv1, fv2 in zip(sc_from_json.fix_node_ids, sc_from_data.fix_node_ids):
        assert_equal(fv1, fv2)
    for vid, spec in sc_from_json.fix_specs.items():
        assert vid in sc_from_data.fix_specs
        assert_equal(spec, sc_from_data.fix_specs[vid])


# @pytest.mark.wip
def test_frame_file_io():
    file_name = 'tower_3D_broken_lines.json'

    here = os.path.dirname(os.path.abspath(__file__))
    file_path = os.path.join(here, 'test_data', file_name)
    node_points, element_vids, fix_node_ids, fix_specs, model_type, material_dict, model_name = \
        read_frame_json(file_path, verbose=True)

    # temp_fp = os.path.join(here, 'tmp.json')
    with TemporaryDirectory() as temp_dir:
        temp_fp = os.path.join(temp_dir, file_name)
        write_frame_json(temp_fp, node_points, element_vids, fix_node_ids, material_dict, 
        fixity_specs=fix_specs, model_type=model_type, model_name = model_name)
        back_node_points, back_element_vids, back_fix_node_ids, back_fix_specs, back_model_type, back_material_dict, back_model_name = \
            read_frame_json(temp_fp, verbose=True)

    for n1, n2 in zip(node_points, back_node_points):
        assert_equal(n1, n2)
    for e1, e2 in zip(element_vids, back_element_vids):
        assert_equal(e1, e2)
    for fv1, fv2 in zip(fix_node_ids, back_fix_node_ids):
        assert_equal(fv1, fv2)
    for vid, spec in fix_specs.items():
        assert vid in back_fix_specs
        assert_equal(spec, back_fix_specs[vid])
    assert model_type == back_model_type
    assert material_dict == back_material_dict
    assert model_name == back_model_name

