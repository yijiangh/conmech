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
    """[summary]
    
    Parameters
    ----------
    json_data : [type]
        [description]
    dim : int, optional
        [description], by default 3
    
    Returns
    -------
    fixities_node_id
        list of int
    fixities_spec
        list of [0-1 specs]
    """
    fixities_spec = {}
    for i, json_node in enumerate(json_data['node_list']):
        if json_node['is_grounded'] == 1:
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
    return fixities_spec


def check_material_dict(mat_dict):
    """check if a material dict is feasible

    example:
    {
    "material_name": "Steel-S235",
    "youngs_modulus": 21000.0,
    "youngs_modulus_unit": "kN/cm2",
    "density": 78.5,
    "density_unit": "kN/m3",
    "poisson_ratio": 0.300149,
    "cross_sec_area": 19.634954084936208,
    "cross_sec_area_unit": "centimeter^2",
    "_cross_sec_area_note": "round_shape, radius 2.5 cm",
    "Jx": 61.35923151542565,
    "Jx_unit": "centimeter^4",
    "_Jx_note": "0.5 * pi * r^4",
    "Iy": 30.679615757712824,
    "Iy_unit": "centimeter^4",
    "_Iy_note": "0.25 * pi * r^4",
    "Iz": 30.679615757712824,
    "Iz_unit": "centimeter^4",
    "_Iz_note": "0.25 * pi * r^4"}

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
            "density" in mat_dict and \
            "density_unit" in mat_dict and \
            "poisson_ratio" in mat_dict and \
            "cross_sec_area" in mat_dict and \
            "Jx" in mat_dict and \
            "Iy" in mat_dict and \
            "Iz" in mat_dict

    # TODO: do other unit conversions
    unit_valid = mat_dict["youngs_modulus_unit"] == 'kN/cm2' and \
                 mat_dict["density_unit"] == 'kN/m3' and \
                 mat_dict["cross_sec_area_unit"] == 'centimeter^2' and \
                 mat_dict["Jx_unit"] == 'centimeter^4' and \
                 mat_dict["Iy_unit"] == 'centimeter^4' and \
                 mat_dict["Iz_unit"] == 'centimeter^4'
    return valid and unit_valid


def extract_model_name_from_path(file_path):
    forename = file_path.split('.json')[0]
    return forename.split(os.sep)[-1]


def read_frame_json(file_path, verbose=False):
    """[summary]
    
    Parameters
    ----------
    file_path : [type]
        [description]
    verbose : bool, optional
        [description], by default False
    
    Returns
    -------
    node_points
        list of [x,y,z]
    element_vids, 
        list of [v_id, v_id]
    fix_node_ids
    
    fix_specs, model_type, material_dicts, model_name
    
    Raises
    ------
    ValueError
        [description]
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

    element_vids = parse_elements(json_data)
    node_points = parse_node_points(json_data, scale=scale, dim=dim)
    fix_specs = parse_fixed_nodes(json_data, dim)

    if 'uniform_cross_section' in json_data and \
       'uniform_material_properties' in json_data and \
       'material_properties' in json_data and \
       json_data["uniform_cross_section"] and \
       json_data["uniform_material_properties"]:
        material_dicts = [json_data['material_properties']] * len(element_vids)
    else:
        material_dicts = [json_element['material_properties'] for json_element in json_data['element_list']]

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
            len(node_points), len(fix_specs), len(element_vids)))
   
    return node_points, element_vids, fix_specs, model_type, material_dicts, model_name, unit


def write_frame_json(file_path, nodes, elements, fixity_specs, material_dicts,
    unif_cross_sec=False, unif_material=False, unit=None, model_type='frame', model_name=None, indent=None, check_material=True):
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
    data['uniform_cross_section'] = unif_cross_sec
    data['uniform_material_properties'] = unif_material
    if unif_cross_sec and unif_material:
        data['material_properties'] = material_dicts[0] if isinstance(material_dicts, list) else material_dicts
        if check_material: assert(check_material_dict(data['material_properties']))
    else:
        data['material_properties'] = {}
        if check_material:
            for mat_dict in material_dicts:
                assert(check_material_dict(mat_dict))
        assert(len(material_dicts) == len(elements))

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
        node_data['is_grounded'] = i in fixity_specs
        if fixity_specs:
            node_data['fixities'] = fixity_specs[i] if node_data['is_grounded'] else []
        else:
            node_data['fixities'] = [1] * 6 if node_data['is_grounded'] else []
        data['node_list'].append(node_data)

    data['element_list'] = []
    for i, element in enumerate(elements):
        element_data = OrderedDict()
        element_data['end_node_ids'] = list([int(v) for v in element])
        element_data['element_id'] = i
        element_data['material_properties'] = {} if unif_cross_sec and unif_material else material_dicts[i]
        data['element_list'].append(element_data)

    with open(file_path, 'w+') as outfile:
        if indent:
            json.dump(data, outfile, indent=indent)
        else:
            json.dump(data, outfile)


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
            assert el_data['description_frame'] == 'global', 'only description in the global frame supported now!'
            uniform_element_load[el_data['applied_element_id']] = [el_data['wx'], el_data['wy'], el_data['wz']]
    
    return point_load, uniform_element_load, include_self_weight