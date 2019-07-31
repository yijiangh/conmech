# Author(s): Andrew Liew (github.com/andrewliew)
# modified by Yijiang Huang (yijiangh@mit.edu)

from __future__ import print_function

import argparse
import os
import json
from pyconmech import stiffness_checker, parse_load_case_from_json

# from compas_fea.cad import rhino
from compas_fea.structure import ElementProperties as Properties
from compas_fea.structure import GeneralStep
from compas_fea.structure import GravityLoad
from compas_fea.structure import PinnedDisplacement
from compas_fea.structure import PointLoad
from compas_fea.structure import Steel
from compas_fea.structure import Structure
from compas_fea.structure import TrussSection

ASSEMBLY_INSTANCE_DIR = os.path.join('..', 'assembly_instances', 'extrusion')
TEST_DATA_DIR = os.path.join('..', 'test_data')

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

    # material property
    p  = 7850
    A1 = 0.0008
    A2 = 0.0005
    A3 = 0.0001

    # Structure
    mdl = Structure(name='truss_frame', path='C:/Temp/')

    # Elements

    # rhino.add_nodes_elements_from_layers(mdl, line_type='TrussElement', layers='elset_main', pL=A1*p)
    # rhino.add_nodes_elements_from_layers(mdl, line_type='TrussElement', layers='elset_diag', pL=A2*p)
    # rhino.add_nodes_elements_from_layers(mdl, line_type='TrussElement', layers='elset_stays', pL=A3*p)
    #
    # # Sets
    #
    # rhino.add_sets_from_layers(mdl, layers=['nset_pins', 'nset_load_v', 'nset_load_h'])
    #
    # # Materials
    #
    # mdl.add(Steel(name='mat_steel', fy=355, p=p))
    #
    # # Sections
    #
    # mdl.add([
    #     TrussSection(name='sec_main', A=A1),
    #     TrussSection(name='sec_diag', A=A2),
    #     TrussSection(name='sec_stays', A=A3),
    # ])
    #
    # # Properties
    #
    # mdl.add([
    #     Properties(name='ep_main', material='mat_steel', section='sec_main', elset='elset_main'),
    #     Properties(name='ep_diag', material='mat_steel', section='sec_diag', elset='elset_diag'),
    #     Properties(name='ep_stays', material='mat_steel', section='sec_stays', elset='elset_stays'),
    # ])
    #
    # # Displacements
    #
    # mdl.add(PinnedDisplacement(name='disp_pinned', nodes='nset_pins'))
    #
    # # Loads
    #
    # mdl.add([
    #     PointLoad(name='load_v', nodes='nset_load_v', z=-15500),
    #     PointLoad(name='load_h', nodes='nset_load_h', x=5000),
    #     GravityLoad(name='load_gravity', elements=['elset_diag', 'elset_main']),
    # ])
    #
    # # Steps
    #
    # mdl.add([
    #     GeneralStep(name='step_bc', displacements=['disp_pinned']),
    #     GeneralStep(name='step_loads', loads=['load_v', 'load_h', 'load_gravity'], factor=1.5, increments=300),
    # ])
    # mdl.steps_order = ['step_bc', 'step_loads']
    #
    # # Summary
    #
    # mdl.summary()
    #
    # # Run
    #
    # mdl.analyse_and_extract(software='abaqus', fields=['u', 's', 'sf', 'cf', 'rf'], ndof=3)
    #
    # rhino.plot_data(mdl, step='step_loads', field='um', radius=0.1, scale=10, cbar_size=0.3)
    # rhino.plot_data(mdl, step='step_loads', field='sxx', radius=0.1, cbar_size=0.3)  # abaqus:sxx opensees:sf1
    # rhino.plot_reaction_forces(mdl, step='step_loads', scale=0.05)
    # rhino.plot_concentrated_forces(mdl, step='step_loads', scale=0.2)
    #
    # print(mdl.get_nodal_results(step='step_loads', field='um', nodes='nset_load_v'))
    # print(mdl.get_nodal_results(step='step_loads', field='rfm', nodes='nset_pins'))
    # print(mdl.get_element_results(step='step_loads', field='sxx', elements='elset_main'))

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
    else:
        file_path = os.path.join(root_directory, ASSEMBLY_INSTANCE_DIR, args.problem)
        load_path = os.path.join(root_directory, ASSEMBLY_INSTANCE_DIR, args.load)

    cm_nD, cm_fR, cm_eR = compute_conmech(file_path, load_path)
    print(cm_nD)
    print(cm_fR)
    print(cm_eR)

if __name__ == '__main__':
    main()
