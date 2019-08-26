"""
Utility functions for reading/writing frame json files.

"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import os
import json
import datetime
from collections import OrderedDict

'''length scale conversion to meter'''
LENGTH_SCALE_CONVERSION = {
    'millimeter': 1e-3,
    'meter': 1.0,
}

def parse_point(json_point, scale=1.0, dim=3):
    if dim == 3:
        return [scale*json_point['X'], scale*json_point['Y'], scale*json_point['Z']]
    else:
        return [scale*json_point['X'], scale*json_point['Y']]


def parse_elements(json_data):
    return [tuple(json_element['end_node_ids']) for json_element in json_data['element_list']]


def parse_node_points(json_data, scale=1.0, dim=3):
    return [parse_point(json_node['point'], scale=scale, dim=dim) for json_node in json_data['node_list']]


def parse_fixed_nodes(json_data, dim=3):
    fixities_node_id = []
    fixities_spec = {}
    for i, json_node in enumerate(json_data['node_list']):
        if json_node['is_grounded'] == 1:
            fixities_node_id.append(i)
            if 'fixities' in json_node and json_node['fixities']:
                if dim == 3:
                    assert len(json_node['fixities']) == dim*2 
                else:
                    assert len(json_node['fixities']) == 3 and dim == 2
                fixities_spec[i] = json_node['fixities']
            else:
                if dim == 3:
                    fixities_spec[i] = [1] * 6
                else:
                    fixities_spec[i] = [1] * 3
    return fixities_node_id, fixities_spec


def check_material_dict(mat_dict):
    """check if a material dict is feasible

    example:
    {
    "material_name": "Steel-S235",
    "youngs_modulus": 21000.0,
    "youngs_modulus_unit": "kN/cm2",
    "shear_modulus": 8076.0,
    "shear_modulus_unit": "kN/cm2",
    "tensile_yeild_stress": 23.5,
    "tensile_yeild_stress_unit": "kN/cm2",
    "density": 78.5,
    "density_unit": "kN/m3",
    "poisson_ratio": 0.300149,
    "radius": 2.5,
    "radius_unit": "centimeter"}

    Parameters
    ----------
    mat_dict : dict
    
    Returns
    -------
    bool
        valid or not
    ------
    """
    valid = "youngs_modulus" in mat_dict and \
            "youngs_modulus_unit" in mat_dict and \
            "shear_modulus" in mat_dict and \
            "shear_modulus_unit" in mat_dict and \
            "density" in mat_dict and \
            "density_unit" in mat_dict and \
            "poisson_ratio" in mat_dict and \
            "radius" in mat_dict and \
            "radius_unit" in mat_dict
    # TODO: element radius should be given to element-level

    # TODO: do other unit conversions
    unit_valid = mat_dict["youngs_modulus_unit"] == 'kN/cm2' and \
                 mat_dict["shear_modulus_unit"] == 'kN/cm2' and \
                 mat_dict["density_unit"] == 'kN/m3' and \
                 mat_dict["radius_unit"] == 'centimeter'
    return valid and unit_valid


def extract_model_name_from_path(file_path):
    forename = file_path.split('.json')[0]
    return forename.split(os.sep)[-1]


def read_frame_json(file_path, verbose=False):
    """Read frame data from a file path.
    
    Parameters
    ----------
    file_path : str
    """
    assert os.path.exists(file_path) and "json file path does not exist!"
    with open(file_path, 'r') as f:
        json_data = json.loads(f.read())
    
    dim = json_data['dimension']
    unit = json_data['unit']
    scale = LENGTH_SCALE_CONVERSION[unit]
    model_type = json_data['model_type']
    model_name = json_data['model_name'] if 'model_name' in json_data else extract_model_name_from_path(file_path)

    if 'frame' in model_type:
        model_type = 'frame'
    elif 'truss' in model_type:
        model_type = 'truss'
    else:
        raise ValueError('model type not supported! Must be frame or truss.')

    material_dict = json_data['material_properties']
    element_vids = parse_elements(json_data)
    node_points = parse_node_points(json_data, scale=scale, dim=dim)
    fix_node_ids, fix_specs = parse_fixed_nodes(json_data, dim)

    if 'node_num' in json_data:
        num_node = json_data['node_num']
        assert len(node_points) == num_node
    if 'element_num' in json_data:
        num_element = json_data['element_num']
        assert len(element_vids) == num_element
    node_id_range = list(range(len(node_points)))
    for ev1, ev2 in element_vids:
        assert ev1 in node_id_range and ev2 in node_id_range, 'element end point id not in node_list id range!'

    if verbose:
        print('Model: {} | Unit: {}'.format(json_data['model_type'], json_data['unit']))
        print('Nodes: {} | Ground: {} | Elements: {}'.format(
            len(node_points), len(fix_node_ids), len(element_vids)))
   
    return node_points, element_vids, fix_node_ids, fix_specs, model_type, material_dict, model_name


def write_frame_json(file_path, nodes, elements, fixed_node_ids, material_dict, 
    fixity_specs={}, unit=None, model_type='frame', model_name=None):
    data = OrderedDict()
    data['model_name'] = model_name if model_name else extract_model_name_from_path(file_path)
    data['model_type'] = model_type
    if not unit:
        print('WARNING: No unit is given in write_frame_json: assuming meter.')
        unit = 'meter'
    else:
        assert unit in LENGTH_SCALE_CONVERSION, 'length unit not supported! please use {}'.format(LENGTH_SCALE_CONVERSION.keys())
    data['unit'] = unit
    data['generate_time'] = str(datetime.datetime.now())
    data['dimension'] = len(nodes[0])
    data['node_num'] = len(nodes)
    data['element_num'] = len(elements)
    assert(check_material_dict(material_dict))
    data['material_properties'] = material_dict

    data['node_list'] = []
    for i, node in enumerate(nodes):
        assert len(node) == data['dimension'], 'node coordinate not in the same dimension!'
        node_data = OrderedDict()
        node_data['point'] = OrderedDict()
        node_data['point']['X'] = node[0] * LENGTH_SCALE_CONVERSION[unit]
        node_data['point']['Y'] = node[1] * LENGTH_SCALE_CONVERSION[unit]
        if data['dimension'] == 3:
            node_data['point']['Z'] = node[2] * LENGTH_SCALE_CONVERSION[unit]
        node_data['node_id'] = i
        node_data['is_grounded'] = i in fixed_node_ids
        if fixity_specs:
            node_data['fixities'] = fixity_specs[i] if node_data['is_grounded'] else []
        else:
            node_data['fixities'] = [1] * 6 if node_data['is_grounded'] else []
        data['node_list'].append(node_data)
    data['element_list'] = []
    for i, element in enumerate(elements):
        element_data = OrderedDict()
        element_data['end_node_ids'] = list(element)
        element_data['element_id'] = i
        data['element_list'].append(element_data)
    with open(file_path, 'w+') as outfile:
        json.dump(data, outfile, indent=4)


def read_load_case_json(file_path):
    """Read load case from a json file.

    Note: 
    - For now, only support `kN` for force unit, `kN-m` for moment unit
    - Gravity is assumed to be along the negative global z axis.
    - Element uniform load is converted to global coordinate in this function
    
    Parameters
    ----------
    file_path : str
    
    Returns
    -------
    point_load : dict 
        {node_id : [Fx, Fy, Fz, Mxx, Myy, Mzz]}, in global coordinate
    uniform_element_load : dict 
        {node_id : [wx, wy, wz]}, in global coordinate
    include_self_weight : bool 
        include self-weight or not, now only supports gravity in global z direction
    """
    assert os.path.exists(file_path) and "json file path does not exist!"
    with open(file_path, 'r') as f:
        json_data = json.loads(f.read())

    # version specific sanity checks
    assert json_data['dimension'] == 3, 'model dimension != 3 not supported now!'
    # TODO: do unit conversion here
    assert json_data['force_unit'] == 'kN'
    assert json_data['moment_unit'] == 'kN-m'

    point_load = {}
    uniform_element_load = {}
    include_self_weight = json_data['include_self_weight'] if 'include_self_weight' in json_data else False
    for pl_data in json_data['point_load_list']:
        point_load[pl_data['applied_node_id']] = [pl_data['Fx'], pl_data['Fy'], pl_data['Fz'], 
                                                 pl_data['Mx'], pl_data['My'], pl_data['Mz']]

    if 'uniformly_distributed_element_load_list' in json_data:
        for el_data in json_data['uniformly_distributed_element_load_list']:
            assert el_data['description_frame'] == 'global'
            uniform_element_load[el_data['applied_element_id']] = [el_data['wx'], el_data['wy'], el_data['wz']]
    
    return point_load, uniform_element_load, include_self_weight

def write_analysis_result_to_json(sc, res_file_path):
    if not sc.has_stored_result():
        print('no result to output!')
        return 

    success, nD, fR, eR = sc.get_solved_results()
    trans_tol, rot_tol = sc.get_nodal_deformation_tol()
    # max_trans, max_rot, max_t_id, max_r_id = sc.get_max_nodal_deformation()
    eR_LG_mats = sc.get_element_local2global_rot_matrices()

    data = OrderedDict()
    data['model_name'] = sc.model_name
    data['model_type'] = sc.model_type
    data['solve_success'] = success
    data['trans_tol'] = trans_tol
    data['rot_tol'] = rot_tol

    data['length_unit'] = 'meter'
    data["force_unit"] =  "kN"
    data["moment_unit"] = "kN-m"
    
    nD_data = []
    for n_id, nd in nD.items():
        nd_data = OrderedDict()
        nd_data['node_id'] = n_id
        nd_data['node_pose'] = list(sc.node_points[n_id])
        nd_data['displacement'] = nd.tolist()
        nD_data.append(nd_data)
    data['node_displacement'] = nD_data

    eR_data = []
    for e_id, er in eR.items():
        er_data = OrderedDict()
        er_data['element_id'] = e_id
        er_data['node_ids'] = sc.elements[e_id]
        er_data['reaction'] = OrderedDict()
        er_data['reaction'][0] = er[0].tolist()
        er_data['reaction'][1] = er[1].tolist()
        eR33 = eR_LG_mats[e_id][:3, :3]
        er_data['local_to_global_transformation'] = eR33.tolist()

        eR_data.append(er_data)
    data['element_reaction'] = eR_data

    fR_data = []
    for v_id, fr in fR.items():
        fr_data = OrderedDict()
        fr_data['node_id'] = v_id
        fr_data['node_pose'] = list(sc.node_points[v_id])
        fr_data['reaction'] = fr.tolist()
        fR_data.append(fr_data)
    data['fixity_reaction'] = fR_data

    with open(res_file_path, 'w+') as outfile:
        json.dump(data, outfile, indent=4)
