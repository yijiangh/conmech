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
    
    unit = json_data['unit']
    # length scale for vertex positions
    scale = LENGTH_SCALE_CONVERSION[unit]
    model_name = json_data['model_name'] if 'model_name' in json_data else extract_model_name_from_path(file_path)

    # * nodal positions
    nodes = parse_nodes(json_data['nodes'], scale=scale)

    # * elements & element tags
    # assume all unspecified elements to be in the tag group ""
    elements = [Element.from_data(e) for e in json_data['elements']]
    element_inds_from_tag = {"":[]}
    for i, e in enumerate(elements):
        if e.elem_tag not in element_inds_from_tag:
            element_inds_from_tag[e.elem_tag] = []
        else:
            element_inds_from_tag[e.elem_tag].append(i)

    # * supports
    supports = [Support.from_data(s) for s in json_data['supports']]

    # unit converted!
    unit = 'meter'

    # * joints
    joints = [Joint.from_data(j) for j in json_data['joints']]
    # elem_tag sanity checks
    for joint in joints:
        if len(joint.elem_tags) == 0:
            joint.elem_tags = [""]
        for e_tag in joint.elem_tags:
            assert e_tag in element_inds_from_tag, 'joint using an element tag not specified in element tag set!'

    # * materials
    materials = [Material.from_data(m) for m in json_data['materials']]
    # elem_tag sanity checks
    for m in materials:
        if len(m.elem_tags) == 0:
            m.elem_tags = [""]
        for e_tag in m.elem_tags:
            assert e_tag in element_inds_from_tag, 'material using an element tag not specified in element tag set!'

    # * cross secs
    crosssecs = [CrossSec.from_data(c) for c in json_data['cross_secs']]
    for cs in crosssecs:
        if len(cs.elem_tags) == 0:
            cs['elem_tags'] = [""]
        for e_tag in cs.elem_tags:
            assert e_tag in element_inds_from_tag, 'cross section using an element tag not specified in element tag set!'

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
        print('Model: {} | Unit: {}'.format(model_name, json_data['unit']))
        print('Nodes: {} | Ground: {} | Elements: {}'.format(
            len(nodes), num_grounded_nodes, len(elements)))
   
    return nodes, elements, supports, joints, materials, crosssecs, model_name, unit


def write_frame_json(file_path, nodes, elements, fixity_specs, material_dicts,
    unif_cross_sec=False, unif_material=False, unit=None, model_type='frame', model_name=None, indent=None, check_material=True):
    raise NotImplementedError()

    # data = OrderedDict()
    # data['model_name'] = model_name if model_name else extract_model_name_from_path(file_path)
    # data['model_type'] = model_type
    # if not unit:
    #     print('WARNING: No unit is given in write_frame_json: assuming meter.')
    #     unit = 'meter'
    # else:
    #     assert unit in LENGTH_SCALE_CONVERSION, 'length unit not supported! please use {}'.format(LENGTH_SCALE_CONVERSION.keys())
    # data['unit'] = unit
    # data['generate_time'] = str(datetime.datetime.now())
    # data['dimension'] = len(nodes[0])
    # data['node_num'] = len(nodes)
    # data['element_num'] = len(elements)
    # data['uniform_cross_section'] = unif_cross_sec
    # data['uniform_material_properties'] = unif_material
    # if unif_cross_sec and unif_material:
    #     data['material_properties'] = material_dicts[0] if isinstance(material_dicts, list) else material_dicts
    #     if check_material: assert(check_material_dict(data['material_properties']))
    # else:
    #     data['material_properties'] = {}
    #     if check_material:
    #         for mat_dict in material_dicts:
    #             assert(check_material_dict(mat_dict))
    #     assert(len(material_dicts) == len(elements))

    # data['node_list'] = []
    # for i, node in enumerate(nodes):
    #     assert len(node) == data['dimension'], 'node coordinate not in the same dimension!'
    #     node_data = OrderedDict()
    #     node_data['point'] = OrderedDict()
    #     node_data['point']['X'] = node[0] * LENGTH_SCALE_CONVERSION[unit]
    #     node_data['point']['Y'] = node[1] * LENGTH_SCALE_CONVERSION[unit]
    #     if data['dimension'] == 3:
    #         node_data['point']['Z'] = node[2] * LENGTH_SCALE_CONVERSION[unit]
    #     node_data['node_id'] = i
    #     node_data['is_grounded'] = i in fixity_specs
    #     if fixity_specs:
    #         node_data['fixities'] = fixity_specs[i] if node_data['is_grounded'] else []
    #     else:
    #         node_data['fixities'] = [1] * 6 if node_data['is_grounded'] else []
    #     data['node_list'].append(node_data)

    # data['element_list'] = []
    # for i, element in enumerate(elements):
    #     element_data = OrderedDict()
    #     element_data['end_node_ids'] = list([int(v) for v in element])
    #     element_data['element_id'] = i
    #     element_data['material_properties'] = {} if unif_cross_sec and unif_material else material_dicts[i]
    #     data['element_list'].append(element_data)

    # with open(file_path, 'w+') as outfile:
    #     if indent:
    #         json.dump(data, outfile, indent=indent)
    #     else:
    #         json.dump(data, outfile)


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