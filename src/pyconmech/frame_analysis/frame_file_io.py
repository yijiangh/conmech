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
import warnings

from pyconmech.frame_analysis.io_base import Node, Element, Support, Joint, CrossSec, Material, PointLoad, UniformlyDistLoad, GravityLoad

'''length scale conversion to meter'''
LENGTH_SCALE_CONVERSION = {
    'millimeter': 1e-3,
    'meter': 1.0,
}

def parse_point(json_point, scale=1.0):
    return [scale*json_point[0], scale*json_point[1], scale*json_point[2]]

def parse_nodes(node_data, scale=1.0):
    nodes = []
    for n in node_data:
        n['point'] = parse_point(n['point'], scale=scale)
        nodes.append(Node.from_data(n))
    return nodes

def check_crosssec_dict(crosssec):
    """check if a material dict is feasible

    Parameters
    ----------
    crosssec : dict
    
    Returns
    -------
    bool
        valid or not
    ------
    """
    if not ("A" in crosssec and "Jx" in crosssec and "Iy" in crosssec and "Iz" in crosssec):
        err_msg = 'Invalid cross section!'
        raise RuntimeError(err_msg)
    # TODO unit checks
    # unit_valid = mat_dict["youngs_modulus_unit"] == 'kN/cm2' and \
    #              mat_dict["density_unit"] == 'kN/m3' and \
    #              mat_dict["cross_sec_area_unit"] == 'cm^2' and \
    #              mat_dict["Jx_unit"] == 'cm^4' and \
    #              mat_dict["Iy_unit"] == 'cm^4' and \
    #              mat_dict["Iz_unit"] == 'cm^4'
    # if not unit_valid:
    #     warnings.warn('Material unit not consistent to the standard one! Be careful!')
    return True

def check_material_dict(material):
    if not ("E" in material and "G12" in material and "density" in material):
        err_msg = 'Invalid material!'
        raise RuntimeError(err_msg)
    return True

def extract_model_name_from_path(file_path):
    forename = file_path.split('.json')[0]
    return forename.split(os.sep)[-1]

def read_frame_json(file_path, verbose=False, strict_check=True):
    """parse a frame structure out of a json file
    
    Note: node point coordinates are converted to meter in this function.

    # TODO merge into Model base class .from_json
    
    Parameters
    ----------
    file_path : [type]
        [description]
    verbose : bool, optional
        [description], by default False
    strict_check : bool, optional
        raise if unit/missing attribute detected, by default True
    
    Returns
    -------
    node_points
        list of [x,y,z], note that the unit is converted to meter
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
    model_name = json_data['model_name'] if 'model_name' in json_data else extract_model_name_from_path(file_path)
    return read_frame_data(json_data, model_name=model_name, verbose=verbose, strict_check=strict_check)

def read_frame_data(json_data, model_name=None, verbose=False, strict_check=True):
    unit = json_data['unit']
    # length scale for vertex positions
    scale = LENGTH_SCALE_CONVERSION[unit]

    # * nodal positions
    nodes = parse_nodes(json_data['nodes'], scale=scale)

    # * elements & element tags
    # ! assume all unspecified elements to be in the tag group ""
    elements = [Element.from_data(e) for e in json_data['elements']]
    element_inds_from_tag = {}
    for e in elements:
        if e.elem_tag not in element_inds_from_tag:
            element_inds_from_tag[e.elem_tag] = []
        element_inds_from_tag[e.elem_tag].append(e.elem_ind)

    # * supports
    supports = [Support.from_data(s) for s in json_data['supports']]

    # unit converted!
    unit = 'meter'

    # * joints
    joints = [Joint.from_data(j) for j in json_data['joints']]

    # * materials
    materials = [Material.from_data(m) for m in json_data['materials']]

    # * cross secs
    crosssecs = [CrossSec.from_data(c) for c in json_data['cross_secs']]

    # * sanity checks
    if 'node_num' in json_data:
        assert len(nodes) == json_data['node_num']
    if 'element_num' in json_data:
        assert len(elements) == json_data['element_num']
    node_id_range = list(range(len(nodes)))
    for e in elements:
        assert e.end_node_inds[0] in node_id_range and \
               e.end_node_inds[1] in node_id_range, \
               'element end point id not in node_list id range!'

    num_grounded_nodes = 0
    for n in nodes:
        if n.is_grounded:
            num_grounded_nodes += 1
    # assert grounded_nodes > 0, 'The structure must have at lease one grounded node!'

    if verbose:
        print('Model: {} | Original unit: {} | Generated time: {}'.format(model_name, json_data['unit'], 
            json_data['generate_time'] if 'generate_time' in json_data else ''))
        print('Nodes: {} | Elements: {} | Supports: {} | Joints: {} | Materials: {} | Cross Secs: {} | Tag Ground: {} '.format(
            len(nodes), len(elements), len(supports), len(joints), len(materials), len(crosssecs), num_grounded_nodes))
   
    return nodes, elements, supports, joints, materials, crosssecs, model_name, unit

def frame_to_data(nodes, elements, supports, joints, materials, crosssecs, model_name):
    data = OrderedDict()
    return data

def read_load_case_json(file_path):
    """Read load case from a json file.

    Note: 
    - For now, only support `kN` for force unit, `kN-m` for moment unit
    - Element uniform load is converted to global coordinate in this function
    
    Parameters
    ----------
    file_path : str
    
    Returns
    -------
    point_load : list 
    uniform_element_load : list 
    gravity_load : object or None 
    """
    assert os.path.exists(file_path) and "json file path does not exist!"
    with open(file_path, 'r') as f:
        json_data = json.loads(f.read())

    point_loads = [PointLoad.from_data(pl) for pl in json_data['ploads']]
    uniform_element_loads = [UniformlyDistLoad.from_data(el) for el in json_data['eloads']]
    gravity_load = None if 'gravity' not in json_data else json_data['gravity']
    
    return point_loads, uniform_element_loads, gravity_load