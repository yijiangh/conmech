from __future__ import print_function
import os
import sys
if sys.version_info[0] < 3:
    from backports import tempfile
else:
    import tempfile

import pytest
import numpy as np
from numpy.linalg import norm, multi_dot
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
    n_attempts=50, existing_ids=[], expect_partial_ids_success=True, trans_tol=2e-3):
    sc = StiffnessChecker(json_file_path=frame_file_path)
    sc.set_loads(point_loads=parsed_pt_loads, include_self_weight=include_sw)
    assert not sc.has_stored_result()
    sc.set_nodal_displacement_tol(trans_tol=trans_tol)

    sol_success = sc.solve(existing_ids, eid_sanity_check=True)

    trans_tol, rot_tol = sc.get_nodal_deformation_tol()
    max_trans, max_rot, max_t_id, max_r_id = sc.get_max_nodal_deformation()
    print('sol benchmark: max_trans {0} / tol {4}, max_rot {1} / tol {5}, max_trans_id {2}, max_rot_id {3}'.format(\
        max_trans, max_rot, max_t_id, max_r_id, trans_tol, rot_tol))

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
        sc_new.set_nodal_displacement_tol(trans_tol=trans_tol)

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
def test_stiffness_checker_solve_consistency(test_case, load_case):
    if test_case == 'tower':
        file_name = 'tower_3D.json'
        load_file_name = 'tower_3D_load_case.json'
        
        # some parts are floating in this partial ids
        failed_existing_ids = [0, 4, 7, 8, 9] 
        if load_case == 'self_weight' or load_case == 'self_weight+point_load':
            success_existing_ids = [1, 12, 5, 10, 4, 11, 13, 3, 2, 0, 7]
            trans_tol = 0.04
        else:
            success_existing_ids = [0,1,2,3,8,9,10,11,12,13,14,15]
            trans_tol = 0.001
    elif test_case == 'topopt-100':
        file_name = 'topopt-100.json'
        load_file_name = 'topopt-100_load_case.json'

        # a big cantilever are floating in this partial ids
        failed_existing_ids = [125, 126, 115, 122, 111, 108, 23, 22, 98, 75, 64, 34, 61, 65, 59, 60, 39, 36, 44, 67]

        success_existing_ids = [0, 1, 2, 3, 4, 5, 6, 8, 11, 14, 15, 19, 20, 21, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33,
                34, 35, 36, 37, 38, 41, 42, 43, 44, 46, 51, 52, 53, 54, 55, 76, 88, 89, 90, 91, 92,
                93, 94, 96, 97, 98, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 
                114, 115, 117, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131]

        trans_tol = 0.002
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
    n_attempts = 100
    
    print('################\nfull solve success checks')
    repetitive_test_stiffness_checker(json_path, parsed_pt_loads=parsed_pt_loads, include_sw=include_sw, \
        n_attempts=n_attempts, existing_ids=[], expect_partial_ids_success=True, trans_tol=trans_tol)

    print('################\npartial solve success checks')
    repetitive_test_stiffness_checker(json_path, parsed_pt_loads=parsed_pt_loads, include_sw=include_sw, \
        n_attempts=n_attempts, existing_ids=success_existing_ids, expect_partial_ids_success=True, trans_tol=trans_tol)

    print('################\npartial solve failure checks')
    repetitive_test_stiffness_checker(json_path, parsed_pt_loads=parsed_pt_loads, include_sw=include_sw, \
        n_attempts=n_attempts, existing_ids=failed_existing_ids, expect_partial_ids_success=False, trans_tol=trans_tol)


def test_init_stiffness_checker():
    file_name = 'tower_3D_broken_lines.json'

    here = os.path.dirname(os.path.abspath(__file__))
    file_path = os.path.join(here, 'test_data', file_name)

    print('init StiffnessChecker from json file path...')
    # this is the same as: sc_from_json = StiffnessChecker(file_path)
    sc_from_json = StiffnessChecker.from_json(file_path, verbose=False)

    print('init StiffnessChecker from frame data...')
    nodes, elements, fixed_node_ids, fix_specs, model_type, material_dicts, model_name = \
        read_frame_json(file_path, verbose=True)
    sc_from_data = StiffnessChecker.from_frame_data(nodes, elements, fixed_node_ids, material_dicts, 
        fixity_specs=fix_specs, unit='meter', model_type=model_type, model_name=model_name, verbose=False)

    assert sc_from_json.model_name == sc_from_data.model_name
    assert sc_from_json.model_type == sc_from_data.model_type
    for mat1, mat2 in zip(sc_from_json.materials, sc_from_json.materials):
        assert mat1 == mat2
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
    node_points, element_vids, fix_node_ids, fix_specs, model_type, material_dicts, model_name = \
        read_frame_json(file_path, verbose=True)

    # temp_fp = os.path.join(here, 'tmp.json')
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_fp = os.path.join(temp_dir, file_name)
        write_frame_json(temp_fp, node_points, element_vids, fix_node_ids, material_dicts, 
        fixity_specs=fix_specs, model_type=model_type, model_name = model_name)
        back_node_points, back_element_vids, back_fix_node_ids, back_fix_specs, back_model_type, back_material_dicts, back_model_name = \
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
    for mat1, mat2 in zip(material_dicts, back_material_dicts):
        assert mat1 == mat2
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


def repetitive_test_equilibrium(frame_file_path, parsed_pt_loads=None, include_sw=False,
        existing_ids=[], expect_partial_ids_success=True):
    debug = False
    sc = StiffnessChecker(json_file_path=frame_file_path)
    sc.set_loads(point_loads=parsed_pt_loads, include_self_weight=include_sw)
    assert not sc.has_stored_result()

    sol_success = sc.solve(existing_ids, eid_sanity_check=True)

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
    for nd in nD.values():
        assert len(nd) == 6
    for fr in fR.values():
        assert len(fr) == 6
    for er in eR.values():
        assert len(er) == 2
        assert len(er[0]) == 6
        assert len(er[1]) == 6
    for fnid in sc.get_element_connected_node_ids(existing_ids=existing_ids, fix_node_only=True):
        assert_equal(nD[fnid], np.zeros(6))

    eR_LG_mats = sc.get_element_local2global_rot_matrices()
    for R in eR_LG_mats.values():
        assert R.shape == (12,12)

    nodal_loads = sc.get_nodal_loads(existing_ids=existing_ids)

    existing_e_ids = existing_ids or list(range(len(sc.elements)))
    existing_node_ids = sc.get_element_connected_node_ids(existing_ids=existing_ids)
    assert len(nodal_loads) == len(existing_node_ids)
    all_node_e_neighbors = sc.get_node_neighbors(return_e_id=True)

    #########################
    # element level checks
    #########################

    print('testing per-element beam equation in local coordinate.')
    # see docs: 
    # https://conmech.readthedocs.io/en/latest/under_the_hood/stiffness_checker_intro.html
    local_stiffness_mats = sc.get_element_stiffness_matrices(in_local_coordinate=True)
    for e_id in existing_e_ids:
        e = sc.elements[e_id]
        K_eL = local_stiffness_mats[e_id]
        eD_G = np.hstack([ nD[e[0]],    nD[e[1]] ])
        eD_L = eR_LG_mats[e_id].dot(eD_G)
        eR_L = np.hstack([ eR[e_id][0], eR[e_id][1] ])
        assert_almost_equal(K_eL.dot(eD_L), eR_L)

    print('testing per-element beam equation in global coordinate.')
    global_stiffness_mats = sc.get_element_stiffness_matrices()
    for e_id in existing_e_ids:
        e = sc.elements[e_id]
        K_eG = global_stiffness_mats[e_id]
        eD_G = np.hstack([ nD[e[0]],    nD[e[1]] ])
        eR_L = np.hstack([ eR[e_id][0], eR[e_id][1] ])
        eR_G = eR_LG_mats[e_id].transpose().dot(eR_L)
        assert_almost_equal(K_eG.dot(eD_G), eR_G)

    #########################
    # node level checks
    #########################

    # test nodal equilibrium, from result elemental reaction force eR
    print('testing nodal equilibrium, from result elemental reaction force')
    for n_id in existing_node_ids:
        internal_force = 0
        if n_id in sc.fix_node_ids:
            fixity_reaction = fR[n_id]
        else:
            fixity_reaction = 0

        connected_e_id_set = all_node_e_neighbors[n_id]
        connected_e_ids = list(connected_e_id_set.intersection(existing_e_ids))
        assert connected_e_ids, 'existing node should at least has one element connected!'
        for e_id in connected_e_ids:
            e = sc.elements[e_id]
            local_id = e.index(n_id) 
            eR_LG = eR_LG_mats[e_id]
            eR_L = np.hstack([ eR[e_id][0], eR[e_id][1] ])
            internal_force += eR_LG.transpose().dot(eR_L)[local_id*6 : local_id*6 + 6]
        
        if debug:
            try:
                assert_almost_equal(fixity_reaction + nodal_loads[n_id], internal_force)
            except AssertionError as err_msg:
                print('----\nnode #{}, fixed: {}, connected elements {},\nfix + load : {},\ninternal_force : {}, \neq_diff {}'.format(
                    n_id, n_id in sc.fix_node_ids, connected_e_ids, 
                    fixity_reaction + nodal_loads[n_id], internal_force,
                    fixity_reaction + nodal_loads[n_id] - internal_force))
                print('\x1b[6;30;43m' + '{}'.format(err_msg) + '\x1b[0m')
        else:
            assert_almost_equal(fixity_reaction + nodal_loads[n_id], internal_force)

    # test nodal force equilibrium, from nodal deformation
    print('testing nodal equilibrium, from nodal deformation')
    local_stiffness_mats = sc.get_element_stiffness_matrices(in_local_coordinate=True)
    global_stiffness_mats = sc.get_element_stiffness_matrices()
    for n_id in existing_node_ids:
        if n_id in sc.fix_node_ids:
            fixity_reaction = fR[n_id]
        else:
            fixity_reaction = 0

        connected_e_id_set = all_node_e_neighbors[n_id]
        connected_e_ids = list(connected_e_id_set.intersection(existing_e_ids))
        e_reaction_node_sum = 0
        for e_id in connected_e_ids:
            e = sc.elements[e_id]
            local_id = e.index(n_id) 

            eD_G = np.hstack([ nD[e[0]], nD[e[1]] ])

            eR_LG = eR_LG_mats[e_id]
            K_eL = local_stiffness_mats[e_id]
            e_reaction = multi_dot([eR_LG.transpose(), K_eL, eR_LG, eD_G])

            K_eG = global_stiffness_mats[e_id]
            e_reaction_rep = multi_dot([K_eG, eD_G])
            assert_almost_equal(e_reaction, e_reaction_rep)
            
            e_reaction_node = e_reaction[local_id*6 : local_id*6 + 6]
            e_reaction_node_sum += e_reaction_node

        if debug:
            try:
                assert_almost_equal(fixity_reaction + nodal_loads[n_id], e_reaction_node_sum)
            except AssertionError as err_msg:
                print('----\nnode #{}, fixed: {}, connected elements {},\nfix + load : {},\ninternal_force : {}, \neq_diff {}'.format(
                    n_id, n_id in sc.fix_node_ids, connected_e_ids, 
                    fixity_reaction + nodal_loads[n_id], e_reaction_node_sum,
                    fixity_reaction + nodal_loads[n_id] - e_reaction_node_sum))
                print('\x1b[6;30;43m' + '{}'.format(err_msg) + '\x1b[0m')
        else:
            assert_almost_equal(fixity_reaction + nodal_loads[n_id], e_reaction_node_sum)


@pytest.mark.equil_check
@pytest.mark.parametrize("test_case, load_case", 
    [('tower', 'self_weight'), ('tower', 'point_load'), ('tower', 'self_weight+point_load'),
     ('topopt-100', 'self_weight'), ('topopt-100', 'point_load'), ('topopt-100', 'self_weight+point_load')])
def test_nodal_equilibrium(test_case, load_case):
    if test_case == 'tower':
        file_name = 'tower_3D.json'
        load_file_name = 'tower_3D_load_case.json'
        
        # some parts are floating in this partial ids
        failed_existing_ids = [0, 4, 7, 8, 9] 
        if load_case == 'self_weight' or load_case == 'self_weight+point_load':
            success_existing_ids = [0, 1, 2, 3, 4, 5, 6, 10, 13]
        else:
            success_existing_ids = [0,1,2,3,8,9,10,11,12,13,14,15]
    elif test_case == 'topopt-100':
        file_name = 'topopt-100.json'
        load_file_name = 'topopt-100_load_case.json'

        # a big cantilever are floating in this partial ids
        failed_existing_ids = [125, 126, 115, 122, 111, 108, 23, 22, 98, 75, 64, 34, 61, 65, 59, 60, 39, 36, 44, 67]
        success_existing_ids = [0, 1, 2, 3, 4, 5, 8, 11, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, \
        43, 44, 51, 52, 53, 88, 89, 90, 91, 92, 93, 94, 96, 97, 98, 102, 103, 104, 105, 106, 107, \
        108, 109, 110, 111, 112, 113, 114, 115, 117, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131]
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
    
    print('\n################\nfull solve success node equilibirum checks')
    repetitive_test_equilibrium(json_path, parsed_pt_loads=parsed_pt_loads, include_sw=include_sw, \
        existing_ids=[], expect_partial_ids_success=True)

    print('\n################\npartial solve success node equilibirum checks')
    repetitive_test_equilibrium(json_path, parsed_pt_loads=parsed_pt_loads, include_sw=include_sw, \
        existing_ids=success_existing_ids, expect_partial_ids_success=True)

    print('\n################\npartial solve failure node equilibirum checks')
    repetitive_test_equilibrium(json_path, parsed_pt_loads=parsed_pt_loads, include_sw=include_sw, \
        existing_ids=failed_existing_ids, expect_partial_ids_success=False)

def repetitive_check_gravity_validity(frame_file_path, existing_ids=[], expect_partial_ids_success=True):
    sc = StiffnessChecker(json_file_path=frame_file_path)

    print('\n################\n [0,0,-100] gravity direction')
    sc.set_loads(include_self_weight=True, gravity_direction=[0,0,-100])
    sc.solve(existing_ids, eid_sanity_check=True)
    compliance_big = sc.get_compliance()

    print('\n################\n [0,0,-1] gravity direction')
    sc.set_loads(include_self_weight=True, gravity_direction=[0,0,-1])
    sc.solve(existing_ids, eid_sanity_check=True)
    compliance_small = sc.get_compliance()

    print('big comp: {}, small comp: {}'.format(compliance_big, compliance_small))
    assert_almost_equal(compliance_big, compliance_small * 1e4)

    success, nD, fR, eR = sc.get_solved_results()

    nodal_loads = sc.get_nodal_loads(existing_ids=existing_ids)
    sw_loads = sc.get_self_weight_loads(existing_ids=existing_ids)
    assert_equal(nodal_loads, sw_loads)

    existing_e_ids = existing_ids or list(range(len(sc.elements)))
    existing_node_ids = sc.get_element_connected_node_ids(existing_ids=existing_ids)
    assert len(nodal_loads) == len(existing_node_ids)
    # all_node_e_neighbors = sc.get_node_neighbors(return_e_id=True)

    node_points, element_vids, _, _, _, material_dicts, _ = \
        read_frame_json(frame_file_path)

    # direct total gravity calculation
    total_weight_force = 0
    for e_id in existing_e_ids:
        n1, n2 = sc.elements[e_id]
        e_len = norm(np.array(node_points[n1]) - np.array(node_points[n2]))

        assert material_dicts[e_id]['cross_sec_area_unit'] == 'centimeter^2' and material_dicts[e_id]["density_unit"] == "kN/m3"
        rho = material_dicts[e_id]['density'] 
        A = material_dicts[e_id]['cross_sec_area'] * 1e-4 # convert to meter^2
        total_weight_force += e_len * A * rho
    
    # total gravity from nodal loads
    total_nodal_weight_load = 0
    for n_id in existing_node_ids:
        assert_equal(nodal_loads[n_id], sw_loads[n_id])
        assert_equal(nodal_loads[n_id][:2], [0, 0])
        total_nodal_weight_load += nodal_loads[n_id][2]
    
    total_z_fixity_reaction = 0
    for n_id, n_fr in fR.items():
        if n_id in existing_node_ids:
            total_z_fixity_reaction += n_fr[2]

    assert_almost_equal(total_weight_force, total_z_fixity_reaction)
    assert_almost_equal(-total_nodal_weight_load, total_z_fixity_reaction)

@pytest.mark.gravity_check
@pytest.mark.parametrize("test_case", 
    [('tower'), ('topopt-100')])
def test_self_weight_validity(test_case):
    if test_case == 'tower':
        file_name = 'tower_3D.json'
        success_existing_ids = list(range(8))
    elif test_case == 'topopt-100':
        file_name = 'topopt-100.json'
        success_existing_ids = [0, 1, 2, 3, 4, 5, 8, 11, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, \
        43, 44, 51, 52, 53, 88, 89, 90, 91, 92, 93, 94, 96, 97, 98, 102, 103, 104, 105, 106, 107, \
        108, 109, 110, 111, 112, 113, 114, 115, 117, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131]
    else:
        assert False, 'not supported test case!'

    root_dir = os.path.dirname(os.path.abspath(__file__))
    json_path = os.path.join(root_dir, 'test_data', file_name)

    print('\n################\nfull solve success gravity validity checks')
    repetitive_check_gravity_validity(json_path, existing_ids=[], expect_partial_ids_success=True)

    print('\n################\npartial solve success gravity validity checks')
    repetitive_check_gravity_validity(json_path, existing_ids=success_existing_ids, expect_partial_ids_success=True)


@pytest.mark.uniform_load_check
@pytest.mark.parametrize("test_case, existing_e_ids", 
    [('tower', []), ('tower', list(range(8))), 
     ('topopt-100', []), ('topopt-100', [0, 1, 2, 3, 4, 5, 8, 11, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, \
        43, 44, 51, 52, 53, 88, 89, 90, 91, 92, 93, 94, 96, 97, 98, 102, 103, 104, 105, 106, 107, \
        108, 109, 110, 111, 112, 113, 114, 115, 117, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131])
     ])
def test_uniformly_distributed_load_with_gravity(test_case, existing_e_ids):
    if test_case == 'tower':
        file_name = 'tower_3D.json'
    elif test_case == 'topopt-100':
        file_name = 'topopt-100.json'
    else:
        assert False, 'not supported test case!'

    root_dir = os.path.dirname(os.path.abspath(__file__))
    json_path = os.path.join(root_dir, 'test_data', file_name)
    sc = StiffnessChecker(json_file_path=json_path)

    # gravity load agreement check
    # _, uniform_element_load, _ = read_load_case_json(load_json_path)
    node_points, element_vids, _, _, _, material_dicts, _ = \
        read_frame_json(json_path)

    uniform_distributed_load = {}
    for e_id in range(len(element_vids)):
        assert material_dicts[e_id]['cross_sec_area_unit'] == 'centimeter^2' and material_dicts[e_id]["density_unit"] == "kN/m3"
        A = material_dicts[e_id]['cross_sec_area'] * 1e-4 # convert to meter^2
        rho = material_dicts[e_id]['density']
        uniform_distributed_load[e_id] = [0, 0, - rho * A]

    # existing_e_ids = list(range(len(sc.elements)))
    sc.set_loads(uniform_distributed_load=uniform_distributed_load, include_self_weight=False)

    nodal_loads = sc.get_nodal_loads(existing_ids=existing_e_ids)
    sw_loads = sc.get_self_weight_loads(existing_ids=existing_e_ids)
    for n_id in nodal_loads.keys():
        assert_almost_equal(nodal_loads[n_id], sw_loads[n_id])

    sc.solve(exist_element_ids=existing_e_ids)
    _, ul_nD, ul_fR, ul_eR = sc.get_solved_results()

    sw_sc = StiffnessChecker(json_file_path=json_path)
    sw_sc.set_loads(include_self_weight=True)
    sw_sc.solve(exist_element_ids=existing_e_ids)
    _, sw_nD, sw_fR, sw_eR = sw_sc.get_solved_results()

    for n_id in ul_nD.keys():
        assert_almost_equal(ul_nD[n_id], sw_nD[n_id])
    for f_id in ul_fR.keys():
        assert_almost_equal(ul_fR[f_id], sw_fR[f_id])
    for e_id in ul_eR.keys():
        assert_almost_equal(ul_eR[e_id][0], sw_eR[e_id][0])
        assert_almost_equal(ul_eR[e_id][1], sw_eR[e_id][1])

# @pytest.mark.ignore(reason='not fully developed')
@pytest.mark.uniform_load_check
@pytest.mark.analy_compare
def test_uniformly_distributed_load_with_analytical_solution():
    """ analytical example in 
            Matrix Structural Analysis 2rd edition, McGuire, Gallagher, Ziemian
            Example 5.8, Page 116 (p.137 in the PDF)
        Material / cross sectional properties from Example 4.8, Page 79 (p.100 in the PDF)

        For more info on the theory of uniformly distributed load lumping, see: 
            Page 108, section 5.2 Loads between nodal points
    """
    file_name = 'uniform_load_analytical_model.json'
    load_file_name = 'uniform_load_analytical_model_load_case.json'

    root_dir = os.path.dirname(os.path.abspath(__file__))
    json_path = os.path.join(root_dir, 'test_data', file_name)
    load_json_path = os.path.join(root_dir, 'test_data', load_file_name)

    sc = StiffnessChecker(json_file_path=json_path)
    point_load, uniform_element_load, include_sw = read_load_case_json(load_json_path)
    sc.set_loads(point_loads=point_load, uniform_distributed_load=uniform_element_load, include_self_weight=include_sw)
    sc.set_nodal_displacement_tol(trans_tol=0.024, rot_tol=0.006)

    def compare_analytical_sol(pass_criteria, nD, fR, eR, nodal_loads, check_decimal=1):
        assert_equal(nodal_loads[1][2], -12.5)
        assert_equal(nodal_loads[1][3], 0.0)
        assert_equal(nodal_loads[1][4], -6.250)
        print('{} \t {}'.format(fR[2][0:3], fR[2][3:6]))

        assert pass_criteria
        assert_equal(nD[0], [0] * 6)
        assert_almost_equal(nD[1], [0, 0, -0.02237, 4.195*1e-3, 5.931*1e-3, 0], decimal=check_decimal)
        assert_equal(nD[2], [0] * 6)
        assert_almost_equal(fR[0], [0, 0, 14.74, -6.45*1e-3, -36.21, 0], decimal=check_decimal)
        # assert_almost_equal(fR[2], [0, 0, 5.25, -41.94, -17.11*1e-3, 0], decimal=check_decimal)

    test_iters = 100
    print('compare analytical res: w/o reinit')
    for _ in range(test_iters):
        sc.solve()
        pass_criteria, nD, fR, eR = sc.get_solved_results()
        nodal_loads = sc.get_nodal_loads()    
        compare_analytical_sol(pass_criteria, nD, fR, eR, nodal_loads)

    print('compare analytical res: w reinit')
    for _ in range(test_iters):
        re_sc = StiffnessChecker(json_file_path=json_path)
        point_load, uniform_element_load, include_sw = read_load_case_json(load_json_path)
        re_sc.set_loads(point_loads=point_load, uniform_distributed_load=uniform_element_load, include_self_weight=include_sw)
        re_sc.set_nodal_displacement_tol(trans_tol=0.024, rot_tol=0.006)

        re_sc.solve()
        pass_criteria, nD, fR, eR = re_sc.get_solved_results()
        nodal_loads = re_sc.get_nodal_loads()    
        compare_analytical_sol(pass_criteria, nD, fR, eR, nodal_loads)


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

