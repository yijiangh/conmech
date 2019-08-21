from __future__ import print_function
import os
from tempfile import TemporaryDirectory

import pytest
import numpy as np
import random
import time
from numpy.testing import assert_equal, assert_almost_equal

from pyconmech import StiffnessChecker
from pyconmech.frame_analysis import read_frame_json, write_frame_json, read_load_case_json, check_material_dict
from pyconmech.database import MATERIAL_PROPERTY_DATABASE


def test_get_material_from_database():
    PLA_MAT = MATERIAL_PROPERTY_DATABASE['PLA']
    STEEL_S235_MAT = MATERIAL_PROPERTY_DATABASE['Steel-S235']
    assert check_material_dict(PLA_MAT)
    assert check_material_dict(STEEL_S235_MAT)


def repetitive_test_stiffness_checker(frame_file_path, parsed_pt_loads=None, include_sw=False,
    n_attempts=50, existing_ids=[], expect_partial_ids_success=True):
    sc = StiffnessChecker(json_file_path=frame_file_path)
    sc.set_loads(point_loads=parsed_pt_loads, include_self_weight=include_sw)
    assert not sc.has_stored_result()

    sol_success = sc.solve(existing_ids, eid_sanity_check=True)

    trans_tol, rot_tol = sc.get_nodal_deformation_tol()
    print('sol benchmark: max_trans {0} / tol {4}, max_rot {1} / tol {5}, max_trans_id {2}, max_rot_id {3}'.format(\
        *sc.get_max_nodal_deformation(), trans_tol, rot_tol))

    if expect_partial_ids_success:
        assert sol_success, 'max_trans {}, max_rot {}, max_trans_id {}, max_rot_id {}'.format(*sc.get_max_nodal_deformation())
    else:
        assert not sol_success, 'max_trans {}, max_rot {}, max_trans_id {}, max_rot_id {}'.format(*sc.get_max_nodal_deformation())

    assert sc.has_stored_result()
    # comparision standard solutions
    success, nD, fR, eR = sc.get_solved_results()
    assert success == sol_success
    assert len(nD) == len(sc.get_element_connected_node_ids(existing_ids=existing_ids))
    if existing_ids:
        assert len(eR) == len(existing_ids)
    else:
        assert len(eR) == len(sc.elements)
    assert len(fR) == len(sc.get_element_connected_node_ids(existing_ids=existing_ids, fix_node_only=True))
    for fnid in sc.get_element_connected_node_ids(existing_ids=existing_ids, fix_node_only=True):
        assert_equal(nD[fnid], np.zeros(6))

    # without reinit
    st_time = time.time()
    for _ in range(n_attempts):
        sol_success = sc.solve(existing_ids)
        tmp_success, tmp_nD, tmp_fR, tmp_eR = sc.get_solved_results()

        assert tmp_success == success and sol_success == success, 'max_trans {}, max_rot {}, max_trans_id {}, max_rot_id {}'.format( \
            *sc.get_max_nodal_deformation())
        for fnid in sc.get_element_connected_node_ids(existing_ids=existing_ids, fix_node_only=True):
            assert_equal(tmp_nD[fnid], np.zeros(6))

        for i in nD.keys() : assert_almost_equal(nD[i], tmp_nD[i]) 
        for i in fR.keys() : assert_almost_equal(fR[i], tmp_fR[i]) 
        for i in eR.keys(): 
            assert_almost_equal(eR[i][0], tmp_eR[i][0]) 
            assert_almost_equal(eR[i][1], tmp_eR[i][1]) 
    tot_time = time.time() - st_time
    print('\nwithout reinit: avg time: {} s, total_time: {} s for {} runs'.format(tot_time / n_attempts, tot_time, n_attempts))

    # reinit
    st_time = time.time()
    for _ in range(n_attempts):
        sc_new = StiffnessChecker(json_file_path=frame_file_path)
        sc_new.set_loads(point_loads=parsed_pt_loads, include_self_weight=include_sw)

        sol_success = sc_new.solve(existing_ids)
        tmp_success, tmp_nD, tmp_fR, tmp_eR = sc_new.get_solved_results()
        assert tmp_success == success and sol_success == success, 'max_trans {}, max_rot {}, max_trans_id {}, max_rot_id {}'.format( \
            *sc_new.get_max_nodal_deformation())
        for fnid in sc_new.get_element_connected_node_ids(existing_ids=existing_ids, fix_node_only=True):
            assert_equal(tmp_nD[fnid], np.zeros(6))

        for i in nD.keys() : assert_almost_equal(nD[i], tmp_nD[i]) 
        for i in fR.keys() : assert_almost_equal(fR[i], tmp_fR[i]) 
        for i in eR.keys(): 
            assert_almost_equal(eR[i][0], tmp_eR[i][0]) 
            assert_almost_equal(eR[i][1], tmp_eR[i][1]) 
    tot_time = time.time() - st_time
    print('\nwith reinit: avg time: {} s, total_time: {} s for {} runs'.format(tot_time / n_attempts, tot_time, n_attempts))


# @pytest.mark.skip(reason=None)
@pytest.mark.parametrize("test_case, load_case", 
    [('tower', 'self_weight'), ('tower', 'point_load'), ('tower', 'self_weight+point_load'),
     ('topopt-100', 'self_weight'), ('topopt-100', 'point_load'), ('topopt-100', 'self_weight+point_load')])
def test_stiffness_checker_consistency_self_weight(test_case, load_case):
    if test_case == 'tower':
        file_name = 'tower_3D.json'
        load_file_name = 'tower_3D_load_case.json'
        
        # some parts are floating in this partial ids
        failed_existing_ids = [0, 4, 7, 8, 9] 
        if load_case == 'self_weight' or load_case == 'self_weight+point_load':
            success_existing_ids = list(range(24))
        else:
            success_existing_ids = [0,1,2,3,8,9,10,11,12,13,14,15]
    elif test_case == 'topopt-100':
        file_name = 'topopt-100.json'
        load_file_name = 'topopt-100_load_case.json'

        # a big cantilever are floating in this partial ids
        failed_existing_ids = [125, 126, 115, 122, 111, 108, 23, 22, 98, 75, 64, 34, 61, 65, 59, 60, 39, 36, 44, 67]
        success_existing_ids = list(range(132)) # full structure for now...
    else:
        assert False, 'not supported test case!'

    root_dir = os.path.dirname(os.path.abspath(__file__))
    json_path = os.path.join(root_dir, 'test_data', file_name)

    if load_case == 'self_weight':
        parsed_pt_loads = []
        include_sw = True
    elif load_case == 'point_load' or load_case == 'self_weight+point_load':
        load_json_path = os.path.join(root_dir, 'test_data', load_file_name)
        parsed_pt_loads, _, _ = read_load_case_json(load_json_path)
        include_sw = load_case=='self_weight+point_load'
    else:
        assert False, 'not supported load case!'
    

    # repetitive tests
    n_attempts = 10
    
    print('################\nfull solve success checks')
    repetitive_test_stiffness_checker(json_path, parsed_pt_loads=parsed_pt_loads, include_sw=include_sw, \
        n_attempts=n_attempts, existing_ids=[], expect_partial_ids_success=True)

    print('################\npartial solve success checks')
    repetitive_test_stiffness_checker(json_path, parsed_pt_loads=parsed_pt_loads, include_sw=include_sw, \
        n_attempts=n_attempts, existing_ids=success_existing_ids, expect_partial_ids_success=True)

    print('################\npartial solve failure checks')
    repetitive_test_stiffness_checker(json_path, parsed_pt_loads=parsed_pt_loads, include_sw=include_sw, \
        n_attempts=n_attempts, existing_ids=failed_existing_ids, expect_partial_ids_success=False)


def test_init_stiffness_checker():
    file_name = 'tower_3D_broken_lines.json'

    here = os.path.dirname(os.path.abspath(__file__))
    file_path = os.path.join(here, 'test_data', file_name)

    print('init StiffnessChecker from json file path...')
    # this is the same as: sc_from_json = StiffnessChecker(file_path)
    sc_from_json = StiffnessChecker.from_json(file_path, verbose=False)

    print('init StiffnessChecker from frame data...')
    nodes, elements, fixed_node_ids, fix_specs, model_type, material_dict, model_name = \
        read_frame_json(file_path, verbose=True)
    sc_from_data = StiffnessChecker.from_frame_data(nodes, elements, fixed_node_ids, material_dict, 
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


def test_get_nodal_loads():
    file_name = 'tower_3D.json'
    load_file_name = 'tower_3D_load_case.json'
    root_dir = os.path.dirname(os.path.abspath(__file__))
    json_path = os.path.join(root_dir, 'test_data', file_name)
    load_json_path = os.path.join(root_dir, 'test_data', load_file_name)

    sc = StiffnessChecker(json_file_path=json_path)
    sc.set_self_weight_load(False)
    parsed_pt_loads, _, _ = read_load_case_json(load_json_path)

    sc.set_loads(point_loads=parsed_pt_loads)
    fetched_nodal_loads = sc.get_nodal_loads()

    for node_id in parsed_pt_loads:
        assert node_id in fetched_nodal_loads
        assert_equal(fetched_nodal_loads[node_id], parsed_pt_loads[node_id])


def test_neighnor_query():
    # TODO: not completed
    file_name = 'tower_3D.json'
    root_dir = os.path.dirname(os.path.abspath(__file__))
    json_path = os.path.join(root_dir, 'test_data', file_name)

    sc = StiffnessChecker(json_file_path=json_path)
    assert_equal(sc.fix_element_ids, [0,1,2,3,8,9,10,11,12,13,14,15])
    assert_equal(sc.get_element_connected_node_ids(fix_node_only=True), [0,1,2,3]) 
    assert_equal(sc.get_element_connected_node_ids(existing_ids=[0,3], fix_node_only=True), [1,2]) 


# @pytest.mark.ii
# def test_helpers():
#     file_name = 'topopt-100.json'
#     root_dir = os.path.dirname(os.path.abspath(__file__))
#     json_path = os.path.join(root_dir, 'test_data', file_name)
#     sc = StiffnessChecker(json_file_path=json_path)
#     ids = []

#     failed_existing_ids = [125, 126, 115, 122, 111, 108, 23, 22, 98, 75, 64, 34, 61, 65, 59, 60, 39, 36, 44, 67]
#     failed_existing_node_ids = sc.get_element_connected_node_ids(failed_existing_ids)
#     print(failed_existing_node_ids)

#     for i, v in enumerate(sc.node_points):
#         if abs(v[2] - 0.155) < 1e-3 and i in failed_existing_node_ids:
#             ids.append(i)
#     print(ids)

