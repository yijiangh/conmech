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
    n_attempts=50, existing_ids=[], expect_partial_ids_success=True):
    sc = StiffnessChecker(json_file_path=frame_file_path)
    sc.set_loads(point_loads=parsed_pt_loads, include_self_weight=include_sw)
    assert not sc.has_stored_result()

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
    with tempfile.TemporaryDirectory() as temp_dir:
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

    from random import shuffle
    # shuffle(existing_e_ids)
    # shuffle(existing_node_ids)

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
        # print('---\nelement #{}: {}, local beam equation diff {}'.format(e_id, e, norm(K_eL.dot(eD_L) - eR_L)))
        # element's internal reaction = elastic deformation caused by nodal displacements (in local coordinate)
        assert_almost_equal(K_eL.dot(eD_L), eR_L)

    print('testing per-element beam equation in global coordinate.')
    global_stiffness_mats = sc.get_element_stiffness_matrices()
    for e_id in existing_e_ids:
        e = sc.elements[e_id]
        K_eG = global_stiffness_mats[e_id]
        eD_G = np.hstack([ nD[e[0]],    nD[e[1]] ])
        eR_L = np.hstack([ eR[e_id][0], eR[e_id][1] ])
        eR_G = eR_LG_mats[e_id].transpose().dot(eR_L)
        # print('---\nelement #{}: {}, global beam equation diff {}'.format(e_id, e, norm(K_eG.dot(eD_G) - eR_G)))
        # element's internal reaction = elastic deformation caused by nodal displacements (in global coordinate)
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
            # continue

        connected_e_id_set = all_node_e_neighbors[n_id]
        connected_e_ids = list(connected_e_id_set.intersection(existing_e_ids))
        assert connected_e_ids, 'existing node should at least has one element connected!'
        # print('before: {}, after: {}'.format(connected_e_id_set, connected_e_ids))
        for e_id in connected_e_ids:
            e = sc.elements[e_id]
            local_id = e.index(n_id) 
            eR_LG = eR_LG_mats[e_id]
            eR_L = np.hstack([ eR[e_id][0], eR[e_id][1] ])
            internal_force += eR_LG.transpose().dot(eR_L)[local_id*6 : local_id*6 + 6]

        # print('node #{}, fixed: {}, connected elements {}, fixity reaction {}, internal force {}, nodal_load {}, eq_diff {}'.format(
        #     n_id, n_id in sc.fix_node_ids, e, 
        #     fixity_reaction, internal_force, nodal_loads[n_id], norm(fixity_reaction + nodal_loads[n_id] - internal_force)))
        
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

    # tower error, eq_diff:
    # node 2, e 0, 13, 14
    # [ 4.44089210e-16 -4.44089210e-16  8.00905782e-01 -2.66968594e-01 -2.66968594e-01 -1.04083409e-17]
    # node 6, e [0, 4, 7, 12, 15, 17, 22]
    # [ 1.13797860e-15  7.21644966e-16  8.00905782e-01  2.66968594e-01 2.66968594e-01 -1.11022302e-16]

    # topopt error eq_diff:
    # node 42, e 0, 89, 88, not fixed
    # [-7.62329653e-21 -1.27054942e-20  2.59944330e-06  6.61744490e-24 -6.61744490e-24 -9.60537572e-25]
    # [-4.33680869e-19  8.67361738e-19  4.62403169e-03 -6.77626358e-21 1.15196481e-19 -2.24195101e-20]
    # node 0, e 0, 131, 91, 93, 94, not fixed
    # [ 1.38966343e-21 -5.08219768e-21  2.59944330e-06 -9.92616735e-24 -2.64697796e-23  8.27180613e-25]
    # [-2.25514052e-17 -1.25767452e-17  4.62403169e-03 -6.77626358e-20 -5.48877350e-19 -4.06575815e-20]

    # test nodal force equilibrium, from nodal deformation
    print('testing nodal equilibrium, from nodal deformation')
    local_stiffness_mats = sc.get_element_stiffness_matrices(in_local_coordinate=True)
    global_stiffness_mats = sc.get_element_stiffness_matrices()
    for n_id in existing_node_ids:
        if n_id in sc.fix_node_ids:
            fixity_reaction = fR[n_id]
        else:
            fixity_reaction = 0
            # continue

        connected_e_id_set = all_node_e_neighbors[n_id]
        connected_e_ids = list(connected_e_id_set.intersection(existing_e_ids))
        # print('before: {}, after: {}'.format(connected_e_id_set, connected_e_ids))
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

        # print('node #{}, fixed: {}, connected elements {}, fixity reaction {}, internal force {}, nodal_load {}, eq_diff {}'.format(
            # n_id, n_id in sc.fix_node_ids, e, 
            # fixity_reaction, internal_force, nodal_loads[n_id], norm(fixity_reaction + nodal_loads[n_id] - internal_force)))

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
# @pytest.mark.parametrize("test_case, load_case", 
#     [('tower', 'self_weight'), ('tower', 'point_load'), ('tower', 'self_weight+point_load'),
#      ('topopt-100', 'self_weight'), ('topopt-100', 'point_load'), ('topopt-100', 'self_weight+point_load')])
@pytest.mark.parametrize("test_case, load_case", 
    # [('topopt-100', 'self_weight'), ('topopt-100', 'point_load'), ('topopt-100', 'self_weight+point_load')])
    [('tower', 'self_weight'), ('topopt-100', 'self_weight')])
def test_nodal_equilibrium(test_case, load_case):
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
    sc.set_loads(include_self_weight=True)

    sol_success = sc.solve(existing_ids, eid_sanity_check=True)

    success, nD, fR, eR = sc.get_solved_results()

    nodal_loads = sc.get_nodal_loads(existing_ids=existing_ids)
    sw_loads = sc.get_self_weight_loads(existing_ids=existing_ids)
    assert_equal(nodal_loads, sw_loads)

    existing_e_ids = existing_ids or list(range(len(sc.elements)))
    existing_node_ids = sc.get_element_connected_node_ids(existing_ids=existing_ids)
    assert len(nodal_loads) == len(existing_node_ids)
    # all_node_e_neighbors = sc.get_node_neighbors(return_e_id=True)

    node_points, element_vids, _, _, _, material_dict, _ = \
        read_frame_json(frame_file_path)
    assert material_dict['radius_unit'] == 'centimeter' and material_dict["density_unit"] == "kN/m3"
    r = material_dict['radius'] * 1e-2 # convert to meter
    rho = material_dict['density'] 
    total_weight_force = 0

    # direct total gravity calculation
    for e_id in existing_e_ids:
        n1, n2 = sc.elements[e_id]
        e_len = norm(np.array(node_points[n1]) - np.array(node_points[n2]))
        # L * A * \rho
        total_weight_force += e_len * np.pi * (r**2) * rho
    
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

    assert_almost_equal(total_weight_force,       total_z_fixity_reaction)
    assert_almost_equal(-total_nodal_weight_load, total_z_fixity_reaction)

@pytest.mark.gravity_check
@pytest.mark.parametrize("test_case", 
    [('tower'), ('topopt-100')])
def test_self_weight_validity(test_case):
    if test_case == 'tower':
        file_name = 'tower_3D.json'
        load_file_name = 'tower_3D_load_case.json'
        
        success_existing_ids = list(range(24))
    elif test_case == 'topopt-100':
        file_name = 'topopt-100.json'
        load_file_name = 'topopt-100_load_case.json'

        # a big cantilever are floating in this partial ids
        success_existing_ids = list(range(132)) # full structure for now...
    else:
        assert False, 'not supported test case!'

    root_dir = os.path.dirname(os.path.abspath(__file__))
    json_path = os.path.join(root_dir, 'test_data', file_name)

    print('\n################\nfull solve success gravity validity hecks')
    repetitive_check_gravity_validity(json_path, existing_ids=[], expect_partial_ids_success=True)

    print('\n################\npartial solve success gravity validity checks')
    repetitive_check_gravity_validity(json_path, existing_ids=success_existing_ids, expect_partial_ids_success=True)


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

