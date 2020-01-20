from pyconmech.frame_analysis import StiffnessChecker

def stiffness_checker_parse_and_solve(file_path, existing_ids=[], pt_loads=None, include_sw=False, result_file_path=None):
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

    if result_file_path:
        sc.write_result_to_json(result_file_path)
    else:
        return sc.result_data