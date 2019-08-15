# Author(s): Andrew Liew (github.com/andrewliew)
# modified by Yijiang Huang (yijiangh@mit.edu)

from __future__ import print_function

import argparse
import os
import json

import numpy as np
from numpy import divide, maximum, abs
from numpy.testing import assert_equal, assert_almost_equal
from pyconmech import stiffness_checker 

from compas_fea.structure import ElementProperties as Properties
from compas_fea.structure import GeneralStep
from compas_fea.structure import GravityLoad
from compas_fea.structure import PointLoad
from compas_fea.structure import ElasticIsotropic
from compas_fea.structure import Structure
from compas_fea.structure import CircularSection
from compas_fea.structure import PinnedDisplacement, GeneralDisplacement

import matplotlib.pyplot as plt

OPENSEES_PATH = 'C:/OpenSees/OpenSees.exe'

ASSEMBLY_INSTANCE_DIR = os.path.join('..', 'assembly_instances', 'extrusion')
TEST_DATA_DIR = os.path.join('..', 'test_data')

__all__ = [
           'run',
           ]

def parse_length_unit_conversion(unit_str):
    """ get the length unit conversion scale to meter (required by compas_fea)

    Parameters
    ----------
    unit_str : str

    Returns
    -------
    scale: double

    """
    if unit_str == 'millimeter':
        scale = 1e-3
    elif unit_str == 'centimeter':
        scale = 1e-2
    elif unit_str == 'meter':
        scale = 1.0
    else:
        raise NotImplementedError('length unit not right!')
    return scale


def parse_pressure_scale_conversion(unit_str):
    """ get the pressure unit conversion to Pascal (required by compas_fea)

    Parameters
    ----------
    unit_str : str

    Returns
    -------
    scale: double

    """
    if unit_str == 'kN/cm2':
        scale = 1e7
    else:
        raise NotImplementedError('pressure unit not right!')
    return scale


def parse_density_scale_conversion(unit_str):
    """ get the pressure unit conversion to kg / m^3 (required by compas_fea)

    Parameters
    ----------
    unit_str : str

    Returns
    -------
    scale: double

    """
    if unit_str == 'kN/m3':
        scale = 101.971621
    else:
        raise NotImplementedError('density unit not right!')
    return scale


def parse_force_scale_conversion(unit_str):
    """ get the force unit conversion to Netwon(N) (required by compas_fea)

    Parameters
    ----------
    unit_str : str

    Returns
    -------
    scale: double

    """
    if unit_str == 'kN':
        scale = 1e3
    else:
        raise NotImplementedError('force unit not right!')
    return scale


def parse_moment_scale_conversion(unit_str):
    """ get the force unit conversion to Netwon x m (N-m) (required by compas_fea)

    Parameters
    ----------
    unit_str : str

    Returns
    -------
    scale: double

    """
    if unit_str == 'kN-m':
        scale = 1e3
    else:
        raise NotImplementedError('moment unit not right!')
    return scale


def parse_point(json_point, scale=1.0):
    return [scale*json_point['X'], scale*json_point['Y'], scale*json_point['Z']]


def parse_frame_nodes(json_data):
    """ parse a list of np array [x, y, z] representing nodal coordinates

    Parameters
    ----------
    json_data : dict

    Returns
    -------
    data: a list of 3-list
        a list of nodal coordinates, unit in meter

    """
    scale = parse_length_unit_conversion(json_data['unit'])
    return [parse_point(json_node['point'], scale)
            for json_node in json_data['node_list']]


def parse_elements(json_data):
    """ parse a list of end node ids for elements

    Parameters
    ----------
    json_data : dict

    Returns
    -------
    data: a list of (int, int) tuple
        a list of nodal id pair, referring to the node list

    """
    return [tuple(json_element['end_node_ids'])
        for json_element in json_data['element_list']]


def parse_fixties(json_data):
    """ parse a list of fixity nodes

    Parameters
    ----------
    json_data : dict

    Returns
    -------
    data: a list of [node_id, fix_x, fix_y, fix_z, fix_mx, fix_my, fix_mz]
        where fix_<> is a bool, =1 means the dof is fixed, =0 means unconstrained.

    """
    data = []
    for i, json_node in enumerate(json_data['node_list']):
        if json_node['is_grounded'] == 1:
            if 'fixities' in json_node:
                fix_dofs = json_node['fixities'] if json_node['fixities'] else [1] * 6
            else:
                fix_dofs = [1] * 6
            data.append([i] + fix_dofs)
    return data


def parse_circular_cross_sec_radius(json_data):
    scale = parse_length_unit_conversion(json_data['material_properties']['radius_unit'])
    return scale *  json_data['material_properties']['radius']



def parse_load_case(json_data):
    """ parse load cases (pt load & gravity)

    Parameters
    ----------
    json_data : str

    Returns
    -------
    pt_loads: list of (7) double list
        each row is [node_Id, Fx, Fy, Fz, Mx, My, Mz], described in global coordiate

    include_sw: bool

    """
    assert 3 == json_data['dimension']
    force_scale = parse_force_scale_conversion(json_data['force_unit'])
    moment_scale = parse_moment_scale_conversion(json_data['moment_unit'])
    pt_loads = []
    include_sw = json_data['include_self_weight']
    for pl in json_data['point_load_list']:
        pt_loads.append([pl['applied_node_id'],
                         force_scale*pl['Fx'],
                         force_scale*pl['Fy'],
                         force_scale*pl['Fz'],
                         moment_scale*pl['Mx'],
                         moment_scale*pl['My'],
                         moment_scale*pl['Mz']])
    return pt_loads, include_sw


def parse_abaqus_result_json(file_name, temp_dir, step='step_loads'):
    res_file_path = os.path.join(temp_dir, file_name, file_name + '-results.json')
    with open(res_file_path, 'r') as f:
        json_data = json.loads(f.read())
        print('abaqus result parsed: {}'.format(res_file_path))

    nodal_data = json_data[step]['nodal']
    n_num = len(nodal_data['ux'])
    assert n_num == len(nodal_data['uy']) and \
        n_num == len(nodal_data['uz']) and \
        n_num == len(nodal_data['urx']) and \
        n_num == len(nodal_data['ury']) and \
        n_num == len(nodal_data['urz'])

    nD = {}
    for n_id in nodal_data['ux'].keys():
        nD[int(n_id)] = [nodal_data['ux'][n_id],
                         nodal_data['uy'][n_id],
                         nodal_data['uz'][n_id],
                         nodal_data['urx'][n_id],
                         nodal_data['ury'][n_id],
                         nodal_data['urz'][n_id]]

    fR = {}
    for n_id in nodal_data['rfx'].keys():
        fR[int(n_id)] = [nodal_data['rfx'][n_id],
                         nodal_data['rfy'][n_id],
                         nodal_data['rfz'][n_id],
                         nodal_data['rmx'][n_id],
                         nodal_data['rmy'][n_id],
                         nodal_data['rmz'][n_id]]

    eR = {}
    return nD, fR, eR


def parse_conmech_result_json(file_name, temp_dir):
    res_file_path = os.path.join(temp_dir, file_name + '-results.json')
    with open(res_file_path, 'r') as f:
        json_data = json.loads(f.read())
        print('conmech result parsed: {}'.format(res_file_path))

    nD_data = json_data['node_displacement']
    nD = {}
    for nd in nD_data:
        # meter / rad, all good
        nD[nd['node_id']] = nd['displacement']

    fR_data = json_data['fixity_reaction']
    fR = {}
    for fr in fR_data:
        fR[fr['node_id']] = np.array(fr['reaction']) * 1e3

    eR = {}
    return nD, fR, eR


def compute_compas_fea(file_path, load_path, fea_engine='abaqus', recompute=True):
    """ Use abaqus (via compas_fea) to perform elastic FEA on the given frame
    under a given load case. If no load path is specified, elemental gravity
    will be assumbed to be applied.

    Parameters
    ----------
    file_path : string
        full path to the frame shape's json file.
    load_path : type
        full path to the load case's json file.

    Returns
    -------
    nD: dict
        Reactional nodal displacement
        key is the node id.
        value is
        (nodal_id, dx, dy, dz, theta_x, theta_y, theta_z).

    fR: dict
        Fixities reaction force, moment.
        key is the nodal id.
        value is [Fxyz, Mxyz] in the global axes.

    eR: dict
        Element-wise reaction force, moment (two ends).
        key is the element id.
        (Fxyz_1, Mxyz_1, Fxyz_2, Mxyz_2)

    """
    root_dir = os.path.dirname(os.path.abspath(__file__))
    temp_dir = os.path.join(root_dir, 'compas_fea-temp')
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    file_json_name = file_path.split(os.sep)[-1]
    file_name = file_json_name.split('.')[0]
    print('compas_fea initing: file name {}'.format(file_name))
    if not recompute:
        nD, fR, eR = parse_abaqus_result_json(file_name, temp_dir)
        return nD, fR, eR

    with open(file_path, 'r') as f:
        json_data = json.loads(f.read())
    load_json_data = {}
    if load_path:
        with open(load_path, 'r') as f:
            load_json_data = json.loads(f.read())

    # init an empty structure
    mdl = Structure(name=file_name, path=os.path.join(temp_dir, ''))

    # nodes
    mdl.add_nodes(nodes=parse_frame_nodes(json_data))

    # elements
    elements = parse_elements(json_data)

    # align local axes with conmech
    sc = stiffness_checker(json_file_path=file_path, verbose=False)
    e_rot_mats = sc.get_element_local2global_rot_matrices()
    assert len(e_rot_mats) == len(elements)
    for e, mat in zip(elements, e_rot_mats):
        # compas_fea local axis convention is differrent to the one used in conmech:
        # in compas_fea
        # 'ex' axis represents the cross-section’s major axis
        # 'ey' is the cross-section’s minor axis
        # 'ez' is the axis along the element
        # TODO: this numpy array to list conversion
        # is essential to make compas_fea work...
        ez = list(mat[0][0:3]) # conmech longitude axis
        ex = list(mat[1][0:3]) # conmech cross sec major axis
        ey = list(mat[2][0:3]) # conmech cross sec minor axis

        mdl.add_element(nodes=e,
                        type='BeamElement',
                        axes={'ex': ex, 'ey': ey, 'ez': ez})
        # print(mdl.elements[mdl.check_element_exists(nodes=e)])

    assert_equal(mdl.element_count(), len(elements))

    # Sets
    # just convenient aliases for referring to a group of elements
    mdl.add_set(name='elset_all', type='element', selection=list(range(mdl.element_count())))

    mdl.add_set(name='nset_all', type='node', selection=list(range(mdl.node_count())))

    fixities = parse_fixties(json_data)
    mdl.add_set(name='nset_fix', type='node', selection=[f[0] for f in fixities])

    if load_json_data:
        pt_loads, include_sw = parse_load_case(load_json_data)
        # mdl.add_set(name='nset_pt_load', type='node', selection=[l[0] for l in pt_loads])
    else:
        pt_loads = []
        include_sw = True
    if pt_loads:
        mdl.add_set(name='nset_v_load_all', type='node', selection=[pl[0] for pl in pt_loads])

    # Materials
    # Young’s modulus E [in units of Pa]
    # Poisson’s ratio v and density p [kg per cubic metre].
    mat_json = json_data['material_properties']
    mat_name = 'mat_' + mat_json['material_name']
    E_scale = parse_pressure_scale_conversion(mat_json['youngs_modulus_unit'])
    p_scale = parse_density_scale_conversion(mat_json['density_unit'])
    mdl.add(ElasticIsotropic(name=mat_name,
                             E=E_scale * mat_json['youngs_modulus'],
                             v=mat_json['poisson_ratio'],
                             p=p_scale * mat_json['density']))

    # G_scale = parse_pressure_scale_conversion(mat_json['shear_modulus_unit'])
    # print('{}, {}'.format(mdl.materials['mat_' + mat_json['material_name']].G, G_scale * mat_json['shear_modulus']))
    # assert_almost_equal(mdl.materials['mat_' + mat_json['material_name']].G['G'], G_scale * mat_json['shear_modulus'])
    # print('-----------material')
    # print(mdl.materials[mat_name])

    # Sections
    # SI units should be used, this includes the use of metres m for cross-section dimensions, not millimetres mm.
    sec_name = 'sec_circ'
    mdl.add(CircularSection(name=sec_name, r=parse_circular_cross_sec_radius(json_data)))

    # print('-----------cross section')
    # print(mdl.sections[sec_name])

    # Properties, associate material & cross sec w/ element sets
    mdl.add(Properties(name='ep_all', material=mat_name, section=sec_name, elset='elset_all'))

    # Displacements
    # pin supports
    for i, fix in enumerate(fixities):
        f_dof = []
        for j in range(6):
            if fix[j+1] == 1:
                f_dof.append(0)
            else:
                f_dof.append(None)
        mdl.add(GeneralDisplacement(name='disp_fix_'+str(i), nodes=[fix[0]], x=f_dof[0], y=f_dof[1], z=f_dof[2], xx=f_dof[3], yy=f_dof[4], zz=f_dof[5]))
    # print('-----------fixities')
    # for i in range(len(fixities)):
    #     print(mdl.displacements['disp_fix_'+str(i)])

    # Loads
    if pt_loads:
        mdl.add([PointLoad(name='load_v_'+str(i), nodes=[pl[0]],
                           x=pl[1], y=pl[2], z=pl[3],
                           xx=pl[4], yy=pl[5], zz=pl[6])
                 for i, pl in enumerate(pt_loads)])
        if include_sw:
            mdl.add(GravityLoad(name='load_gravity', elements='elset_all'))
    else:
        mdl.add(GravityLoad(name='load_gravity', elements='elset_all'))
    # print('-----------loads')
    # print(mdl.loads['load_gravity'])
    # for i in range(len(pt_loads)):
    #     print(mdl.loads['load_v_'+str(i)])

    # Steps
    loads_names = []
    if pt_loads:
        loads_names.extend(['load_v_'+str(i) for i in range(len(pt_loads))])
    if include_sw:
        loads_names.append('load_gravity')
    mdl.add([
        GeneralStep(name='step_bc', displacements=['disp_fix_'+str(i) for i in range(len(fixities))]),
        GeneralStep(name='step_loads', loads=loads_names)
        ])

    # a boundary condition step such as 'step_bc' above, should always be applied as the first step to prevent rigid body motion
    mdl.steps_order = ['step_bc', 'step_loads']

    # Summary
    mdl.summary()

    # Run
    # node
    # 'u': nodal displacement: ux, uy, uz, um (magnitude)
    # 'ur': nodal rotation
    # 'rf': reaction force
    # 'cf': concentrated force (external load)
    # 'cm': concentrated moment (external load)

    # element
    # 's': beam stress (conmech cannot compute this at
    # version 0.1.1)
    # For beam, the following values are evaluated
    # at the "integration point" 'ip1' (middle point)
    # and pts along the axis: 'sp3, sp7, sp11, sp15'
    # sxx: axial
    # syy: hoop
    # sxy: torsion
    # smises: Von Mises
    # smaxp: max principal
    # sminp: min principal

    # 'sf': beam section force
    # sf1: axial
    # sf2: shear x
    # sf3: shear y
    if fea_engine == 'abaqus':
        mdl.analyse_and_extract(software='abaqus', fields=['u', 'ur', 'rf', 'rm', 'sf'], ndof=6, output=True)
        nD, fR, eR = parse_abaqus_result_json(file_name, temp_dir)
    elif fea_engine == 'opensees':
        mdl.analyse_and_extract(software='opensees', fields=['u'], exe=OPENSEES_PATH, ndof=6, output=True, save=True)
        raise NotImplementedError('opensees from compas_fea is not fully supported at this moment...')

        nD = {}
        fR = {}
        eR = {}
        # nD = mdl.get_nodal_results(step='step_load', field='ux', nodes='nset_all')
        print(mdl.results)
    else:
        raise NotImplementedError('FEA engine not supported!')

    return nD, fR, eR

def compute_conmech(file_path, load_path='', recompute=True):
    """ Use pyconmech to perform elastic FEA on the given frame under a given
    load case. If no load path is specified, elemental gravity will be assumbed
    to be applied.

    Parameters
    ----------
    file_path : string
        full path to the frame shape's json file.
    load_path : type
        full path to the load case's json file.

    Returns
    -------
    nodal_disp: (n_element x 7) array
        Reactional nodal displacement. Each row is
        (nodal_id, dx, dy, dz, theta_x, theta_y, theta_z).

    fixities_reaction: (n_fixities x 7)
        Fixities reaction force, moment. The first entry of each row is
        nodal_id.

    element_reaction: (n_element x 13) array
        Element-wise reaction force, moment (two ends). The first entry of each
        row is the element id.

    """
    root_dir = os.path.dirname(os.path.abspath(__file__))
    temp_dir = os.path.join(root_dir, 'conmech-temp')
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    file_json_name = file_path.split(os.sep)[-1]
    file_name = file_json_name.split('.')[0]
    print('conmech initing: file name {}'.format(file_name))
    if not recompute:
        nD, fR, eR = parse_conmech_result_json(file_name, temp_dir)
        return nD, fR, eR

    sc = stiffness_checker(json_file_path=file_path, verbose=False)

    sc.set_output_json(True)
    sc.set_output_json_path(file_path = temp_dir, file_name = file_name + "-results.json")

    if load_path:
        ext_load, include_sw = parse_load_case_from_json(load_path)
        sc.set_load(nodal_forces = ext_load)
        sc.set_self_weight_load(include_sw)
        print('conmech using load from: {}, pt_load #{}, include_sw: {}'.format(load_path, len(ext_load), include_sw))
    else:
        print('No external load specified, default using gravity load.')
        sc.set_self_weight_load(True)

    sc.solve()
    success, nD, fR, eR  = sc.get_solved_results()

    # print('============================')
    # print("conmech pass criteria?\n {0}".format(success))
    # trans_tol, rot_tol = sc.get_nodal_deformation_tol()
    # max_trans, max_rot, max_trans_vid, max_rot_vid = sc.get_max_nodal_deformation()
    # compliance = sc.get_compliance()
    #
    # print('max deformation: translation: {0} / tol {1}, at node #{2}'.format(max_trans, trans_tol, max_trans_vid))
    # print('max deformation: rotation: {0} / tol {1}, at node #{2}'.format(max_rot, rot_tol, max_rot_vid))
    # print('compliance: {}'.format(compliance))

    nD, fR, eR = parse_conmech_result_json(file_name, temp_dir)
    return nD, fR, eR


def compute_relative_error(np_vec1, np_vec2, use_abs=True):
    # rel_err = np.true_divide(abs(np_vec1 - np_vec2), maximum(abs(np_vec1), abs(np_vec2)))
    if use_abs:
        rel_err = np.true_divide(abs(np_vec1 - np_vec2), abs(np_vec2))
    else:
        rel_err = np.true_divide(np_vec1 - np_vec2, abs(np_vec2))

    zero_eps = 1e-30
    # if both entries are zero, rel error 0
    for i in range(len(rel_err)):
        if abs(np_vec1[i]) < zero_eps and abs(np_vec2[i]) < zero_eps:
            rel_err[i] = 0.0
    # return np.amax(rel_err)
    return rel_err

def run():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--problem', default='tower_3D.json', help='The name of the problem file to solve')
    parser.add_argument('-fea', '--fea_engine', default='abaqus', help='fea engines name, default to abaqus, avaible: opensees')
    parser.add_argument('-l', '--load', default='tower_3D_load_case.json', help='The name of the load case file to solve')
    parser.add_argument('-a', '--assembly', action='store_true', help='Use test problems from assembly_instance')
    parser.add_argument('-pa', '--parse_result', action='store_true', help='do not recompute abaqus, parse existing results')
    # parser.add_argument('-pl', '--plot', action='store_true', help='plot results')

    args = parser.parse_args()
    print('Arguments:', args)

    root_directory = os.path.dirname(os.path.abspath(__file__))
    file_path = ''
    load_path = ''
    if not args.assembly:
        file_path = os.path.join(root_directory, TEST_DATA_DIR, args.problem)
        if args.load:
            load_path = os.path.join(root_directory, TEST_DATA_DIR, args.load)
    else:
        file_path = os.path.join(root_directory, ASSEMBLY_INSTANCE_DIR, args.problem)

    cm_nD, cm_fR, cm_eR = compute_conmech(file_path, load_path, recompute=not args.parse_result)
    print('===================')
    print('---conmech result:')
    for n_id, nd in cm_nD.items():
        print('node u #{}: {}'.format(n_id, nd))
    print('--------------------')

    # for i in range(len(cm_fR)):
    #     print('fix node #{} r: {}'.format(int(cm_fR[i][0]), cm_fR[i][1:7]))
    # print('--------------------')

    # for i in range(len(cm_fR)):
    #     print('fix node #{} n0: {}'.format(int(cm_fR[i][0]), cm_fR[i][1:7]))
    #     print('element r #{} n1: {}'.format(int(cm_fR[i][0]), cm_fR[i][8:14]))
    #

    ab_nD, ab_fR, ab_eR = compute_compas_fea(file_path, load_path, recompute=not args.parse_result, fea_engine=args.fea_engine)
    # print('===================')
    # print('---abaqus result:')
    # for n_id, nd in ab_nD.items():
    #     print('node u #{}: {}'.format(n_id, nd))
    # print('--------------------')

    # align nodal translational displacement
    cm_nd_trans = []
    ab_nd_trans = []
    cm_nd_rot = []
    ab_nd_rot = []
    for n_id in cm_nD.keys():
        cm_nd_trans.extend(cm_nD[n_id][0:3])
        ab_nd_trans.extend(ab_nD[n_id][0:3])
        cm_nd_rot.extend(cm_nD[n_id][3:6])
        ab_nd_rot.extend(ab_nD[n_id][3:6])

    trans_err = compute_relative_error(np.array(cm_nd_trans), np.array(ab_nd_trans), use_abs=False)
    rot_err = compute_relative_error(np.array(cm_nd_rot), np.array(ab_nd_rot), use_abs=False)
    # print('nodal displ trans rel error: {}, rot rel error: {}'.format(trans_err, rot_err))

    nD_fig, nD_axes = plt.subplots(2, 3)

    ax_name = {0: 'x', 1: 'y', 2: 'z', 3: 'xx', 4: 'yy', 5: 'zz'}
    nDt_len = len(trans_err)

    for i in range(0, 3):
        nD_axes[0, i].plot(list(range(0,len(cm_nD))), trans_err[i:nDt_len:3])
        # nD_axes[0, i].set_xlabel('node id')
        # nD_axes[0, i].set_ylabel('{}: (cm_nD - ab_nD) / ab_nD'.format(ax_name[i]))
        nD_axes[0, i].set_title('nD_{} relative error'.format(ax_name[i]))

        nD_axes[1, i].plot(list(range(0,len(cm_nD))), rot_err[i:nDt_len:3])
        # nD_axes[1, i].set_xlabel('node id')
        # nD_axes[1, i].set_ylabel('{}: (cm_nD - ab_nD) / ab_nD'.format(ax_name[i]))
        nD_axes[1, i].set_title('nD_{} relative error'.format(ax_name[i+3]))

    nD_fig.show()

    # align fixities reaction
    cm_fr_force = []
    ab_fr_force = []
    cm_fr_moment = []
    ab_fr_moment = []
    for n_id in cm_fR.keys():
        cm_fr_force.extend(cm_fR[n_id][0:3])
        ab_fr_force.extend(ab_fR[n_id][0:3])
        cm_fr_moment.extend(cm_fR[n_id][3:6])
        ab_fr_moment.extend(ab_fR[n_id][3:6])
    fR_force_err = compute_relative_error(np.array(cm_fr_force), np.array(ab_fr_force), use_abs=False)
    fR_moment_err = compute_relative_error(np.array(cm_fr_moment), np.array(ab_fr_moment), use_abs=False)
    print('fixities reaction force rel error: {}, moment rel error: {}'.format(fR_force_err, fR_moment_err))

    fR_fig, fR_axes = plt.subplots(2, 3)

    ax_name = {0: 'Fx', 1: 'Fy', 2: 'Fz', 3: 'Mx', 4: 'My', 5: 'Mz'}
    fRt_len = len(fR_force_err)

    for i in range(0, 3):
        fR_axes[0, i].plot(list(range(0,len(cm_fR))), fR_force_err[i:fRt_len:3])
        # fR_axes[0, i].set_xlabel('node id')
        # fR_axes[0, i].set_ylabel('{}: (cm_fR - ab_fR) / ab_fR'.format(ax_name[i]))
        fR_axes[0, i].set_title('fR_{} relative error'.format(ax_name[i]))

        fR_axes[1, i].plot(list(range(0,len(cm_fR))), fR_moment_err[i:fRt_len:3])
        # fR_axes[1, i].set_xlabel('node id')
        # fR_axes[1, i].set_ylabel('{}: (cm_fR - ab_fR) / ab_fR'.format(ax_name[i]))
        fR_axes[1, i].set_title('fR_{} relative error'.format(ax_name[i+3]))

    fR_fig.show()
    print('type <Enter> to exit... ')
    input()


if __name__ == '__main__':
    run()
