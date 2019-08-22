from pyconmech.frame_analysis import StiffnessChecker

def stiffness_checker_parse_file_and_save_result(file_path, result_file_path, existing_ids=[], pt_loads=None, include_sw=False):
    """Wrapper function call to create a StiffnessChecker, solve, and return results
    Mainly for the rpc call, since it doesn't support class creation.
    
    Parameters
    ----------
    file_path : [type]
        [description]
    existing_ids : list, optional
        [description], by default []
    pt_loads : [type], optional
        [description], by default None
    include_sw : bool, optional
        [description], by default False
    
    Returns
    -------
    [type]
        [description]
    """
    sc = StiffnessChecker(json_file_path=file_path)
    sc.set_loads(point_loads=pt_loads, include_self_weight=include_sw)
    sc.solve(existing_ids, eid_sanity_check=True)

    success, nD, fR, eR = sc.get_solved_results()

    trans_tol, rot_tol = sc.get_nodal_deformation_tol()
    # max_trans, max_rot, max_t_id, max_r_id = sc.get_max_nodal_deformation()

    eR_LG_mats = sc.get_element_local2global_rot_matrices()

    sc.write_result_to_json(result_file_path)

    # return (success, nD, fR, eR), (trans_tol, rot_tol), eR_LG_mats