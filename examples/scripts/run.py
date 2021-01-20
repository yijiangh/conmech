import os
import argparse
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
LOAD_EPS = 1e-3

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

def read_supports(file_path, nodes, align_axis=None):
    with open(file_path, 'r') as fp:
        txt_lines = fp.readlines()
        bc_points = []
        bc_cond = []
        for i, line in enumerate(txt_lines):
            txt_sp = line.split(',')
            coord = list(map(float, txt_sp[:3]))
            cond = list(map(int, txt_sp[3:]))
            bc_points.append(coord)
            bc_cond.append(cond)
    supps = []
    if align_axis is None:
        for pt, cond in zip(bc_points, bc_cond):
            node_ind = None
            for i, node in enumerate(nodes):
                if norm(np.array(pt) - np.array(node.point)) < BC_EPS:
                    node_ind = i
                    break
            # else:
            #     raise ValueError('Support point {} not found'.format(pt))
            if node_ind is not None:
                # assume rotational dof is restrained as well at supports
                supps.append(Support(cond + [1,1,1], node_ind))
    else:
        val = np.array(bc_points)[0,align_axis]
        assert all([abs(v-val)<EPS for v in np.array(bc_points)[:,align_axis]]), \
            'Given BC points do not have identical {} axis values'.format(align_axis)
        for i, node in enumerate(nodes):
            if abs(node.point[align_axis]-val) < EPS:
                supps.append(Support([1,1,1] + [1,1,1], i))
    return supps

def read_point_loads(file_path, nodes):
    with open(file_path, 'r') as fp:
        txt_lines = fp.readlines()
        pl_points = []
        pl_vectors = []
        for i, line in enumerate(txt_lines):
            txt_sp = line.split(',')
            coord = list(map(float, txt_sp[:3]))
            load_v = list(map(float, txt_sp[3:]))
            pl_points.append(coord)
            pl_vectors.append(load_v)
    point_loads = []
    for pt, force in zip(pl_points, pl_vectors):
        node_ind = None
        for i, node in enumerate(nodes):
            if norm(np.array(pt) - np.array(node.point)) < LOAD_EPS:
                node_ind = i
                break
        # else:
        #     raise ValueError('Support point {} not found'.format(pt))
        if node_ind is not None:
            moment = [0,0,0]
            point_loads.append(PointLoad(force, moment, node_ind, loadcase=0))
    return point_loads

################################################

def analyze_truss(problem, viewer=True, align_axis=None):
    if align_axis is not None:
        assert align_axis in [0,1,2]
    problem_path = os.path.join(HERE, problem)

    # * parse nodes from txt file
    # ! we assume the unit is meter
    nodes = read_nodes(os.path.join(problem_path, NODE_FILE_NAME))

    # * parse elements from txt file
    elems = read_elems(os.path.join(problem_path, ELEM_FILE_NAME))

    # * parse supports from txt file
    supps = read_supports(os.path.join(problem_path, BC_FILE_NAME), nodes, align_axis=align_axis)

    # * parse loads from txt file
    point_loads = read_point_loads(os.path.join(problem_path, LOAD_FILE_NAME), nodes)

    # * assuming uniform material for all elements
    # Steel, kN/m2
    E = 210000000.0
    G12 = 80760000.0
    # material strength in the specified direction (local x direction)
    fy = 235000.0
    density = 78.5 # kN/m3
    materials = [Material(E, G12, fy, density)]

    # * assuming uniform cross sections for all elements
    A = 0.006 # m2
    Jx = 3E-07 # m4
    Iy = Iz = 0.0002 # m4
    cross_secs = [CrossSec(A, Jx, Iy, Iz)]

    # * assemble info into a model
    joints = []
    model = Model(nodes, elems, supps, joints, materials, cross_secs,  model_name=problem)

    sc = StiffnessChecker(model, checker_engine="numpy", verbose=True)
    loadcase = LoadCase(point_loads=point_loads)
    sc.set_loads(loadcase)
    
    # if the structure's nodal deformation exceeds 1e-3 meter,
    # we want the checker to return `sol_success = False`
    sc.set_nodal_displacement_tol(trans_tol=np.inf, rot_tol=np.inf)
    
    sol_success = sc.solve()
    
    # Get all the analysis information:
    # nodal deformation, fixity reactions, element reactions
    success, nD, fR, eR = sc.get_solved_results()

    # nD is the nodal displacements
    # a dictionary indexed by nodal indices, values are the [ux, uy, uz, theta_x, theta_y, theta_z] deformation vector
    
    # Compliance is the elastic energy: https://en.wikipedia.org/wiki/Elastic_energy
    # the higher this value is, the more flexible a structure is
    # we want this value to be low
    compliance = sc.get_compliance()

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
        # draw supports
        support_pts = np.array([nodes[s.node_ind].point for s in supps])
        ax.scatter3D(support_pts[:,0], support_pts[:,1], support_pts[:,2], c='brown')
        # draw loads
        for pl in point_loads:
            node_point = nodes[pl.node_ind].point
            force = np.array(pl.force)
            force *= 0.05/norm(force)
            ax.quiver(*node_point, *force, lw=1, color='blue')

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
    parser.add_argument('--bc_align_axis', type=int, default=None,
                        help='axis value used to determine BC.')
    args = parser.parse_args()

    analyze_truss(args.problem, viewer=args.viewer, align_axis=args.bc_align_axis)

# Issue:
# python .\examples\scripts\run.py -v --bc_align_axis 0

if __name__ == '__main__':
    main()