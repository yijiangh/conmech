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
    return [json_element['end_node_ids'] for json_element in json_data['element_list']]


def parse_node_points(json_data, scale=1.0, dim=3):
    return [parse_point(json_node['point'], scale=scale, dim=dim) for json_node in json_data['node_list']]


def parse_fixed_nodes(json_data, dim=3):
    fixities_node_id = []
    fixities_spec = {}
    for i, json_node in enumerate(json_data['node_list']):
        if json_node['is_grounded'] == 1:
            fixities_node_id.append(i)
            if json_node['fixities']:
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

