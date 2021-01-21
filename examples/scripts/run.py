import os
import json
import argparse
from termcolor import cprint
import numpy as np
from  numpy.linalg import norm as norm
from pyconmech.frame_analysis import StiffnessChecker
from pyconmech.frame_analysis import Model, Node, Element, Support, Material, CrossSec
from pyconmech.frame_analysis import PointLoad, GravityLoad, UniformlyDistLoad, LoadCase

from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))

NODE_FILE_NAME = 'node.txt'
ELEM_FILE_NAME = 'elem.txt'
LOAD_FILE_NAME = 'load.txt'
BC_FILE_NAME = 'BC.txt'

EPS = 1e-8

# BC point searching radius
BC_EPS = 1e-3
# load point searching radius
LOAD_EPS = 1e-2

################################################

def read_nodes(file_path):
    with open(file_path, 'r') as fp:
        txt_lines = fp.readlines()
        nodes = []
        for i, line in enumerate(txt_lines):
            coord = list(map(float, line.split(',')))
            is_grounded = False # change later
            nodes.append(Node(coord, i, is_grounded))
    return nodes

def read_elems(file_path, index_base=1):
    with open(file_path, 'r') as fp:
        txt_lines = fp.readlines()
        elems = []
        for i, line in enumerate(txt_lines):
            end_node_ids = list(map(lambda x: int(x)-index_base, line.split(',')))
            elems.append(Element(end_node_ids, i, elem_tag='', bending_stiff=True))
    return elems

def read_supports(file_path):
    supps = []
    with open(file_path, 'r') as fp:
        txt_lines = fp.readlines()
        for i, line in enumerate(txt_lines):
            txt_sp = line.split(',')
            cond = list(map(int, txt_sp))
            if not all([c == 0 for c in cond]):
                supps.append(Support(cond + [1,1,1], i))
    return supps

def read_point_loads(file_path):
    point_loads = []
    moment = [0,0,0]
    with open(file_path, 'r') as fp:
        txt_lines = fp.readlines()
        for i, line in enumerate(txt_lines):
            txt_sp = line.split(',')
            load_v = list(map(float, txt_sp))
            if not all([c == 0 for c in load_v]):
                point_loads.append(PointLoad(load_v, moment, i, loadcase=0))
    return point_loads

################################################

def A_solid_cir(radius):
    return np.pi * radius**2

def Jx_solid_cir(radius):
    # https://en.wikipedia.org/wiki/Polar_moment_of_inertia
    return np.pi * radius**4 / 2

def Iy_solid_cir(radius):
    # https://www.engineeringtoolbox.com/area-moment-inertia-d_1328.html
    return np.pi * radius**4 / 4

def solid_cir_crosssec(r, elem_tags=None):
    A = A_solid_cir(r)
    Jx = Jx_solid_cir(r)
    Iy = Iy_solid_cir(r)
    return CrossSec(A, Jx, Iy, Iy, elem_tags=elem_tags, family='Circular', name='solid_circle')

################################################

def analyze_truss(problem, viewer=True, save_model=False, exagg=1.0):
    problem_path = os.path.join(HERE, problem)

    # * parse nodes from txt file
    # ! we assume the unit is meter
    nodes = read_nodes(os.path.join(problem_path, NODE_FILE_NAME))

    # * parse elements from txt file
    elems = read_elems(os.path.join(problem_path, ELEM_FILE_NAME))

    # * parse supports from txt file
    supps = read_supports(os.path.join(problem_path, BC_FILE_NAME))

    # * parse loads from txt file
    point_loads = read_point_loads(os.path.join(problem_path, LOAD_FILE_NAME))

    # * assuming uniform material for all elements
    # Steel, kN/m2
    E = 210000000.0
    G12 = 80760000.0
    # material strength in the specified direction (local x direction)
    fy = 235000.0
    density = 78.5 # kN/m3
    materials = [Material(E, G12, fy, density, family='Steel', name='fake')]

    # * assuming uniform cross sections for all elements
    radius = 0.002 # meter
    cross_secs = [solid_cir_crosssec(radius)]

    # * assemble info into a model
    joints = []
    model = Model(nodes, elems, supps, joints, materials, cross_secs,  model_name=problem)
    loadcase = LoadCase(point_loads=point_loads)
    if save_model:
        save_path = os.path.join(problem_path, problem + '.json')
        with open(save_path, 'w') as f:
            json.dump(model.to_data(), f, indent=None)
        cprint('Model saved to {}'.format(save_path), 'green')

        lc_save_path = os.path.join(problem_path, problem + '_loadcase.json')
        with open(lc_save_path, 'w') as f:
            json.dump(loadcase.to_data(), f, indent=None)
        cprint('Load Case saved to {}'.format(lc_save_path), 'green')

    sc = StiffnessChecker(model, checker_engine="numpy", verbose=True)
    sc.set_loads(loadcase)
    
    # if the structure's nodal deformation exceeds 1e-3 meter,
    # we want the checker to return `sol_success = False`
    sc.set_nodal_displacement_tol(trans_tol=np.inf, rot_tol=np.inf)
    
    success = sc.solve()
    
    # Get all the analysis information:
    # nodal deformation, fixity reactions, element reactions
    success, nD, fR, eR = sc.get_solved_results()
    cprint('Solve success: {}'.format(success), 'green' if success else 'red')

    # nD is the nodal displacements
    # a dictionary indexed by nodal indices, values are the [ux, uy, uz, theta_x, theta_y, theta_z] deformation vector
    trans_tol, rot_tol = sc.get_nodal_deformation_tol()
    max_trans, max_rot, max_trans_vid, max_rot_vid = sc.get_max_nodal_deformation()
    # translation = np.max(np.linalg.norm([d[:3] for d in nD.values()], ord=2, axis=1))
    # rotation = np.max(np.linalg.norm([d[3:] for d in nD.values()], ord=2, axis=1))
    cprint('Max translation deformation: {0:.5f} m / {1:.5} = {2:.5}, at node #{3}'.format(
        max_trans, trans_tol, max_trans / trans_tol, max_trans_vid), 'cyan')
    cprint('Max rotation deformation: {0:.5f} rad / {1:.5} = {2:.5}, at node #{3}'.format(
        max_rot, rot_tol, max_rot / rot_tol, max_rot_vid), 'cyan')
    
    # Compliance is the elastic energy: https://en.wikipedia.org/wiki/Elastic_energy
    # The inverse of stiffness is flexibility or compliance
    # the higher this value is, the more flexible a structure is
    # we want this value to be low
    compliance = sc.get_compliance()
    cprint('Elastic energy: {}'.format(compliance), 'cyan')

    # volume can be computed by simply summing over `cross sectional area * bar length`

    if viewer:
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        # draw elems
        for e in elems:
            n1, n2 = e.end_node_inds
            x = [nodes[n1].point[0], nodes[n2].point[0]]
            y = [nodes[n1].point[1], nodes[n2].point[1]]
            z = [nodes[n1].point[2], nodes[n2].point[2]]
            ax.plot3D(x, y, z, 'black', linewidth=0.5)
        # draw deformed elems
        for e in elems:
            n1, n2 = e.end_node_inds
            e_id = e.elem_ind
            x = np.array([nodes[n1].point[0] + exagg * nD[n1][0], nodes[n2].point[0] + exagg * nD[n2][0]])
            y = np.array([nodes[n1].point[1] + exagg * nD[n1][1], nodes[n2].point[1] + exagg * nD[n2][1]])
            z = np.array([nodes[n1].point[2] + exagg * nD[n1][2], nodes[n2].point[2] + exagg * nD[n2][2]])
            axial_e_reaction = eR[e_id][0][0]
            color = 'blue' if axial_e_reaction < 0 else 'red'
            ax.plot3D(x, y, z, c=color, linewidth=0.5)
        # draw supports
        support_pts = np.array([nodes[s.node_ind].point for s in supps])
        ax.scatter3D(support_pts[:,0], support_pts[:,1], support_pts[:,2], c='brown')
        # draw loads
        for pl in point_loads:
            node_point = nodes[pl.node_ind].point
            force = np.array(pl.force)
            force *= 0.05/norm(force)
            ax.quiver(*node_point, *force, lw=1, color='purple')

        ax.set_title('Truss')
        ax.set_ylabel("Y")
        ax.set_xlabel("X")
        ax.set_zlabel("Z")
        plt.show()

#####################################################

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--problem', default='cantilever_beam_truss',
                        help='The name of the problem to solve')
    parser.add_argument('-v', '--viewer', action='store_true', 
                        help='Enables the viewer, default False')
    parser.add_argument('-s', '--save', action='store_true', 
                        help='Save conmech model, default False')
    parser.add_argument('--exagg', type=float, default=10.0,
                        help='Deformation exaggeration ratio.')
    args = parser.parse_args()

    analyze_truss(args.problem, viewer=args.viewer, save_model=args.save, exagg=args.exagg)

# Issue:
# python .\examples\scripts\run.py -v

if __name__ == '__main__':
    main()