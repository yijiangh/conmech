# Author(s): Andrew Liew (github.com/andrewliew)
# modified by Yijiang Huang (yijiangh@mit.edu)

from __future__ import print_function

import argparse
import os
import json

from numpy.testing import assert_equal, assert_almost_equal
from pyconmech import stiffness_checker, parse_load_case_from_json

from compas_fea.structure import ElementProperties as Properties
from compas_fea.structure import GeneralStep
from compas_fea.structure import GravityLoad
from compas_fea.structure import PointLoad
from compas_fea.structure import ElasticIsotropic
from compas_fea.structure import Structure
from compas_fea.structure import CircularSection
from compas_fea.structure import PinnedDisplacement, GeneralDisplacement

ASSEMBLY_INSTANCE_DIR = os.path.join('..', 'assembly_instances', 'extrusion')
TEST_DATA_DIR = os.path.join('..', 'test_data')

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
            fix_dofs = json_node['fixities'] if json_node['fixities'] else [1] * 6
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


def compute_abaqus(file_path, load_path):
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
    nodal_disp: (n_element x 7) array
        Reactional nodal displacement. Each row is
        (nodal_id, dx, dy, dz, theta_x, theta_y, theta_z).

    element_reaction: (n_element x 13) array
        Element-wise reaction force, moment (two ends). The first entry of each
        row is the element id.

    fixities_reaction: (n_fixities x 7)
        Fixities reaction force, moment. The first entry of each row is
        nodal_id.

    """
    with open(file_path, 'r') as f:
        json_data = json.loads(f.read())
    load_json_data = {}
    if load_path:
        with open(load_path, 'r') as f:
            load_json_data = json.loads(f.read())

    root_dir = os.path.dirname(os.path.abspath(__file__))
    temp_dir = os.path.join(root_dir, 'compas_fea-temp-')
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)
    file_json_name = file_path.split(os.sep)[-1]
    file_name = file_json_name.split('.')[0]
    print('compas_fea initing: file name {}'.format(file_name))

    # init an empty structure
    mdl = Structure(name=file_name, path=temp_dir)

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
    mdl.add(GravityLoad(name='load_gravity', elements='elset_all'))
    if pt_loads:
        mdl.add([PointLoad(name='load_v_'+str(i), nodes=[pl[0]],
                           x=pl[1], y=pl[2], z=pl[3],
                           xx=pl[4], yy=pl[5], zz=pl[6])
                 for i, pl in enumerate(pt_loads)])
    # print('-----------loads')
    # print(mdl.loads['load_gravity'])
    # for i in range(len(pt_loads)):
    #     print(mdl.loads['load_v_'+str(i)])

    # Steps
    mdl.add([
        GeneralStep(name='step_bc', displacements=['disp_fix_'+str(i) for i in range(len(fixities))]),
        GeneralStep(name='step_loads', loads=['load_v_'+str(i) for i in range(len(pt_loads))] + ['load_gravity'])
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
    # 'cf': concentrated force
    # 'cm': concentrated moment

    # element
    # 's': beam stress: sxx, syy, szz
    # 'sf': section force: sf1, sf2, sf3

    mdl.analyse_and_extract(software='abaqus', fields=['u', 'ur', 's', 'sf', 'cf', 'rf'], ndof=6)

    print(mdl.get_nodal_results(step='step_loads', field='um', nodes='nset_v_load_all'))
    print(mdl.get_nodal_results(step='step_loads', field='rfm', nodes='nset_fix'))
    print(mdl.get_element_results(step='step_loads', field='sxx', elements='elset_all'))

    nodal_disp = []
    fixities_reaction = []
    element_reaction = []
    return nodal_disp, fixities_reaction, element_reaction

def compute_conmech(file_path, load_path=''):
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
    sc = stiffness_checker(json_file_path=file_path, verbose=False)
    if load_path:
        ext_load, include_sw = parse_load_case_from_json(load_case_path)
        sc.set_load(nodal_forces = ext_load)
        sc.set_self_weight_load(include_sw)
    else:
        sc.set_self_weight_load(True)

    sc.solve()
    success, nodal_disp, fixities_reaction, element_reaction  = sc.get_solved_results()

    print('============================')
    print("conmech pass criteria?\n {0}".format(success))
    trans_tol, rot_tol = sc.get_nodal_deformation_tol()
    max_trans, max_rot, max_trans_vid, max_rot_vid = sc.get_max_nodal_deformation()
    compliance = sc.get_compliance()

    print('max deformation: translation: {0} / tol {1}, at node #{2}'.format(max_trans, trans_tol, max_trans_vid))
    print('max deformation: rotation: {0} / tol {1}, at node #{2}'.format(max_rot, rot_tol, max_rot_vid))
    print('compliance: {}'.format(compliance))

    return nodal_disp, fixities_reaction, element_reaction


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--problem', default='tower_3D.json', help='The name of the problem file to solve')
    parser.add_argument('-l', '--load', default='tower_3D_load_case.json', help='The name of the load case file to solve')
    parser.add_argument('-a', '--assembly', action='store_true', help='Use test problems from assembly_instance')

    args = parser.parse_args()
    print('Arguments:', args)

    root_directory = os.path.dirname(os.path.abspath(__file__))
    file_path = ''
    load_path = ''
    if not args.assembly:
        file_path = os.path.join(root_directory, TEST_DATA_DIR, args.problem)
        load_path = os.path.join(root_directory, TEST_DATA_DIR, args.load)
    else:
        file_path = os.path.join(root_directory, ASSEMBLY_INSTANCE_DIR, args.problem)

    # cm_nD, cm_fR, cm_eR = compute_conmech(file_path, load_path)
    ab_nD, ab_fR, ab_eR = compute_abaqus(file_path, load_path)


if __name__ == '__main__':
    main()
