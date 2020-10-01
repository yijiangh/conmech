import warnings
import numpy as np
import scipy
import scipy.sparse.linalg as SPLA
import numpy.linalg as LA
from collections import defaultdict
from numpy.linalg import norm
from scipy.sparse import csc_matrix
from pyconmech.frame_analysis.stiffness_base import StiffnessBase
from pyconmech.frame_analysis.io_base import mu2G, G2mu
from termcolor import cprint

DOUBLE_EPS = 1e-14
DOUBLE_INF = 1e14
DEFAULT_GRAVITY_DIRECTION = np.array([0.,0.,-1.])

############################################

# internal calculation unit
LENGTH_UNIT = 'm' # meter
ROT_ANGLE_UNIT = 'rad'
FORCE_UNIT = 'kN'

# derived units
MOMENT_UNIT = FORCE_UNIT + "*" + ROT_ANGLE_UNIT
AREA_UNIT = LENGTH_UNIT + "2"
AREA_INERTIA_UNIT = LENGTH_UNIT + "4"
PRESSURE_UNIT = FORCE_UNIT + "/" + LENGTH_UNIT + "2"; # "kN/m^2" modulus
DENSITY_UNIT  = FORCE_UNIT + "/" + LENGTH_UNIT + "3"; # "kN/m^3"

LENGTH_CONVERSION = {
    'meter' : 1.0,
    'millimeter' : 1e-3,
}

############################################

def axial_stiffness_matrix(L, A, E, cr1=np.inf, cr2=np.inf):
    """local stiffness matrix for a bar with only axial force at two ends

            AE/L (u1 - u2) = F1
            AE/L (u2 - u1) = F2

        Virtual displacement (principle of virtual work) derivation:
        (same result can be derived using variational method or Galerkin method,
        which generalizes FEA to problems beyond structural mechanics.)

        We use linear (C0) interpolation between (x1, phi1) and ((x1, phi2))
        Thus the polynimial basis is [1, x]'
        The interpolation polynomial gives:
            [phi1, phi2]' = A [a0, a1]'
        where A = [[1, x1], [1, x2]]
        Inverting A gives:
            Ainv = (1/(x2-x1)) [[x2, -x1], [-1, 1]]
        Thus the *shape function* N:
            N = [(x2-x)/(x2-x1), (x-x1)/(x2-x1)]

        Thus, the displacement u can be interpolated over an element:
            u = N d
        where d is the nodal displacement dof of the element.

        Strains are determined from displacements:
            \epsilon = B d
        where [B] = [\partial][N] = d/dx [(L-x)/L, x/L]= [-1/L, L]

        Thus, the stiffness matrix 
            [k] = \int_0^L [B]^T E [B] A dx = AE [[1, -1], [-1, 1]]
    
    Parameters
    ----------
    L : [type]
        bar element length
    A : [type]
        cross section area
    E : [type]
        Young's modulus
    cr1: float
        rotational stiffness at the first end of the beam, default to infinite (rigid joint)
    cr2: float
        rotational stiffness at the second end of the beam, default to infinite (rigid joint)

    Return
    ------
    K_ax_x : 2x2 numpy array
    """
    K = np.ones([2,2])
    K[0,1] = -1.
    K[1,0] = -1.
    d_inv = 0.0
    if abs(cr1) < DOUBLE_EPS or abs(cr2) < DOUBLE_EPS:
        # cr1 = cr2 = 0, fully released
        d_inv = 0.0
    else:
        cr1_inv = 1/cr1 if cr1<DOUBLE_INF else 0.0
        cr2_inv = 1/cr2 if cr2<DOUBLE_INF else 0.0
        d_inv = 1.0 / (L/(E*A) + cr1_inv + cr2_inv) 
    K *= d_inv
    return K


def torsional_stiffness_matrix(L, J, G, cr1=np.inf, cr2=np.inf):
    """[summary]

    [Mx1, Mx2]' = ((GJ)/L) * [[1,-1], [-1,1]] [theta_x1, theta_x2]
    
    Parameters
    ----------
    L : [type]
        bar element length
    J : [type]
        torsional constant, unit: length unit^4
        In the case of circular, cylindrical shaft, it's equal to
        the polar moment of inertia of the cross section.
    G : [type]
        modulus of rigidity
    cr1: float
        rotational stiffness at the first end of the beam, default to infinite (rigid joint)
    cr2: float
        rotational stiffness at the second end of the beam, default to infinite (rigid joint)

    Return
    ------
    K_tor_x : 2x2 numpy array
    """
    return axial_stiffness_matrix(L, J, G, cr1=cr1, cr2=cr2)


def bending_stiffness_matrix(L, E, Iz, axis=2, cr1=np.inf, cr2=np.inf):
    """stiffness matrix of uniform beam element without the transverse shear deformation

    Here the stress-strain (\sigma ~ \epsilon) turns into 
    bending moment-curvature (M ~ \kappa).

        strain e_x = - y/\kappa = - y (d^2 v / d^2 x)
    where \kappa is the radius of curvature

    Then from \sigma_x = E e_x,
        \sigma_x = - E y (d^2 v / d^2 x)
    
    From moment equilibrium,
        M_z = - \int_A \sigma_x y dA = E Iz (d^2 v / d^2 x)

    where the lateral displacement v(x) = [N]{d} and \kappa = [B] {d}
        nodal dof {d} = [v_1, theta_z1, v2, theta_z2]
        [B] = d^2/d^2 x [N]

    We use cubic curve interpolation (C1) and [N]'s each row represents
    the Langrange's interpolation functions (polynomials of deg 3 in this case).
    Thus, the stiffness matrix [k] = int_0^L (E Iz d^2[N]) dx

    Parameters
    ----------
    L : [type]
        [description]
    E : [type]
        [description]
    Iz : [type]
        moment of inertia of the section about the z axis
            Iz = \int_A y^2 dA
    axis: int
        1 = local y axis , 2 = local z axis, default to 2
    cr1: float
        rotational stiffness at the first end of the beam, default to infinite (rigid joint)
    cr2: float
        rotational stiffness at the second end of the beam, default to infinite (rigid joint)

    Return
    ------
    K_bend_z : 4x4 numpy array
    """
    K = np.zeros([4,4])
    sign = 1. if axis==2 else -1.
    assert abs(E) > DOUBLE_EPS and abs(Iz) > DOUBLE_EPS
    a1 = cr1*L/E*Iz # k1 = a1 * EI/L
    a2 = cr2*L/E*Iz # k2 = a2 * EI/L
    # a = a1*a2 / (a1*a2 + 4*a1 + 4*a2 + 12)
    # aa = a * (1 + (a1+a2)/(a1*a2))
    # a2_2 = a * (1 + 2/a2)
    # a2_3 = a * (1 + 3/a2)
    # a1_2 = a * (1 + 2/a1)
    # a1_3 = a * (1 + 3/a1)
    if abs(a1) < DOUBLE_INF and abs(a2) < DOUBLE_INF:
        a_demo = a1*a2 + 4*a1 + 4*a2 + 12
        a = a1*a2 / a_demo
        aa = (a1*a2 + a1 +a2) / a_demo
        a2_2 = (a2 + 2)*a1 / a_demo
        a2_3 = (a2 + 3)*a1 / a_demo
        a1_2 = (a1 + 2)*a2 / a_demo
        a1_3 = (a1 + 3)*a2 / a_demo
    elif abs(a1) < DOUBLE_EPS and abs(a2) > DOUBLE_INF:
        # a1=0, a2 = inf
        # a = a1 / (a1 + 4)
        a = 0.0
        aa = 0.25
        a2_2 = 0.0
        a2_3 = 0.0
        a1_2 = 0.5
        a1_3 = 0.75
    elif abs(a2) < DOUBLE_EPS and abs(a1) > DOUBLE_INF:
        a = 0.0
        aa = 0.25
        a1_2 = 0.0
        a1_3 = 0.0
        a2_2 = 0.5
        a2_3 = 0.75
    else:
        assert abs(a1) > DOUBLE_EPS and abs(a2) > DOUBLE_EPS
        inv_a1 = 0.0 if abs(a1) > DOUBLE_INF else 1/a1
        inv_a2 = 0.0 if abs(a2) > DOUBLE_INF else 1/a2
        inv_a_demo = 1+4*inv_a1+4*inv_a2+12*inv_a1*inv_a2
        a = 1 / inv_a_demo
        aa = (1+inv_a1+inv_a2) / inv_a_demo
        a2_2 = (1+2*inv_a1) / inv_a_demo
        a2_3 = (1+3*inv_a1) / inv_a_demo
        a1_2 = (1+2*inv_a2) / inv_a_demo
        a1_3 = (1+3*inv_a2) / inv_a_demo

    # * the original rigid end formula
    # K[0,:] = np.array([12.,      sign*6*L, -12.,      sign*6*L])
    # K[1,:] = np.array([sign*6*L, 4*(L**2), sign*-6*L, 2*(L**2)])
    # K[2,:] = -K[0,:]
    # K[3,:] = np.array([sign*6*L, 2*(L**2), -sign*6*L, 4*(L**2)])
    # K *= (E*Iz/L**3)

    # * the modified formula for joint rotational stiffness
    # See [GMZ] p.394
    K[0,:] = np.array([12./L**2*aa, sign*6/L*a2_2, -12./L**2*aa, sign*6/L*a1_2])
    K[1,:] = np.array([sign*6/L*a2_2, 4*a2_3, sign*(-6/L)*a2_2, 2*a])
    K[2,:] = -K[0,:]
    K[3,:] = np.array([sign*6/L*a1_2, 2*a, sign*(-6/L)*a1_2, 4*a1_3])
    K *= E*Iz/L

    return K

def create_local_stiffness_matrix(L, A, Jx, Iy, Iz, E, mu, ct1=[np.inf]*3, ct2=[np.inf]*3, cr1=[np.inf]*3, cr2=[np.inf]*3):
    """complete 12x12 stiffness matrix for a bisymmetrical member.

        Since for small displacements the axial force effects, torsion,
        and bending about each axis are uncoupled, the influence coeff
        **relating** these effects are zero.
    
    Parameters
    ----------
    L : float
        element length
    A : float
        cross section area
    Jx : float
        torsional constant, unit: length unit^4
        In the case of circular, cylindrical shaft, it's equal to
        the polar moment of inertia of the cross section.
    Iy : float
        moment of inertia \int_A z^2 dA
    Iz : float
        moment of inertia \int_A y^2 dA
    E : float
        Young's modulus
    mu : float
        Poisson ratio
    ct1 : list of three float
        [Tx, Ty, Tz]
        translational spring stiffness at joint 1, set to 0 will release the corresponding dof
        default to np.inf (rigid joint)
    ct2 : list of three float
        [Tx, Ty, Tz]
        translational spring stiffness at joint 2, set to 0 will release the corresponding dof
        default to np.inf (rigid joint)
    cr1 : list of three float
        [Rx, Ry, Rz]
        rotational spring stiffness at joint 1, set to 0 will release the corresponding dof
        default to np.inf (rigid joint)
    cr2 : list of three float
        [Rx, Ry, Rz]
        rotational spring stiffness at joint 2, set to 0 will release the corresponding dof
        default to np.inf (rigid joint)
    
    Returns
    -------
    K : 12x12 numpy array
    """
    assert ct1[1] > DOUBLE_INF and ct2[1] > DOUBLE_INF, 'local y axis translational dof release not implemented!'
    assert ct1[2] > DOUBLE_INF and ct2[2] > DOUBLE_INF, 'local z axis translational dof release not implemented!'
    for i in range(3):
        assert abs(cr1[i]) < DOUBLE_EPS or abs(cr1[i]) > DOUBLE_INF, 'rotational stiffness has some bugs, not fully supported yet'
        assert abs(cr2[i]) < DOUBLE_EPS or abs(cr2[i]) > DOUBLE_INF, 'rotational stiffness has some bugs, not fully supported yet'

    G = mu2G(mu, E)
    # Fx1, Fx2 : u1, u2
    # 0, 6 : 0, 6
    axial_x_k = axial_stiffness_matrix(L, A, E, cr1=ct1[0], cr2=ct2[0])
    # Mx1, Mx2 : \theta_x1, \theta_x2
    # 3, 9 : 3, 9
    tor_x_k = torsional_stiffness_matrix(L, Jx, G, cr1=cr1[0], cr2=cr2[0])

    # Fy1, Mz1, Fy2, Mz2 : v1, \theta_z1, v2, \theta_z2
    # 1, 5, 7, 11 : 1, 5, 7, 11
    bend_z_k = bending_stiffness_matrix(L, E, Iz, axis=2, cr1=cr1[2], cr2=cr2[2])
    # Fz1, My1, Fz2, My2 : v1, \theta_z1, v2, \theta_z2
    # 2, 4, 8, 10 : 2, 4, 8, 10
    bend_y_k = bending_stiffness_matrix(L, E, Iz, axis=1, cr1=cr1[1], cr2=cr2[1])

    K = np.zeros([12,12])
    K[np.ix_([0,6], [0,6])] += axial_x_k
    K[np.ix_([3,9], [3,9])] += tor_x_k
    K[np.ix_([1,5,7,11], [1,5,7,11])] += bend_z_k
    K[np.ix_([2,4,8,10], [2,4,8,10])] += bend_y_k
    return K

def global2local_transf_matrix(end_vert_u, end_vert_v, rot_y2x=0.0):
    assert len(end_vert_u) == len(end_vert_v)
    assert len(end_vert_u) == 2 or len(end_vert_u) == 3
    dim = len(end_vert_u)
    L = norm(end_vert_u-end_vert_v)
    assert L > 1e-6

    # by convention, the new x axis is along the element's direction
    # directional cosine of the new x axis in the global world frame
    c_x = (end_vert_v[0] - end_vert_u[0])/L
    c_y = (end_vert_v[1] - end_vert_u[1])/L
    R = np.zeros([3,3])

    if rot_y2x != 0.0:
        raise NotImplementedError('rotation around local axis not implemented!')

    if 3 == dim:
        c_z = (end_vert_v[2] - end_vert_u[2]) / L
        # TODO rotaxis
        if abs(abs(c_z) - 1.0) < DOUBLE_EPS:
            # the element is parallel to global z axis
            # cross product is not defined, in this case
            # it's just a rotation about the global z axis
            # in x-y plane
            R[0, 2] = -c_z
            R[1, 1] = 1
            R[2, 0] = c_z
        else:
            # local x_axis = element's vector
            new_x = np.array([c_x, c_y, c_z])
            # local y axis = cross product with global z axis
            new_y = -np.cross(new_x, [0,0,1.0])
            new_y /= norm(new_y)
            new_z = np.cross(new_x, new_y)
            R[0, :] = new_x
            R[1, :] = new_y
            R[2, :] = new_z
        # R = R * rot_axis;
    elif 2 == dim:
        R[0,:] = [c_x, c_y, 0]
        R[1,:] = [-c_y, c_x, 0]
        R[2,2] = 1.0
        # TODO check rotaxis
    return R

def turn_diagblock(R3):
    """compile a 4x4 block-diagonal matrix

    Parameters
    ----------
    R3 : 3x3 numpy array

    Returns
    -------
    12x12 numpy array
    """
    R = np.zeros((12,12))
    for i in range(4):
        R[i*3:(i+1)*3, i*3:(i+1)*3] = R3
    return R

###############################################

def assemble_global_stiffness_matrix(k_list, nV, id_map, exist_e_ids=[]):
    """assemble element stiffness matrix into a sparse system stiffness matrix
    
    Parameters
    ----------
    k_list : list of np array
        a list of equal-sized global stiffness matrix
    nV : int
        number of vertices, used to compute total dof for sizing 
        the output sparse matrix
    id_map : np array
        (nEx(2xnode_dof)) idex matrix:
            M[element index] = [nodal dof index]
        e.g. for a single beam
            M[0, :] = [0, 1, 2, 3, 4, 5]

    Return
    ------
    Ksp : scipy.sparse.csc_matrix
        https://docs.scipy.org/doc/scipy-0.15.1/reference/generated/scipy.sparse.csc_matrix.html#scipy.sparse.csc_matrix
    """
    node_dof = 6
    total_dof = node_dof * nV

    if len(exist_e_ids)==0:
        exist_e_ids = range(len(k_list))
    else:
        assert set(exist_e_ids) <= set(range(len(k_list)))

    row = []
    col = []
    data = []
    for e_id in exist_e_ids:
        Ke = k_list[e_id]
        assert Ke.shape == (node_dof*2, node_dof*2)
        for i in range(node_dof*2):
            for j in range(node_dof*2):
                if abs(Ke[i,j]) > DOUBLE_EPS:
                    row.append(id_map[e_id, i])
                    col.append(id_map[e_id, j])
                    data.append(Ke[i,j])
    Ksp = csc_matrix((data, (row, col)), shape=(total_dof, total_dof), dtype=float)
    return Ksp

###############################################

def compute_lumped_uniformly_distributed_load(w_G, R_LG, L):
    """Use fixed-end beam formula to lump uniformly distributed load 
    to equivalent nodal loads
    
    Parameters
    ----------
    w_G : np array of size (3,)
        uniformly distributed load vector in the global coordinate
        unit: force unit / length unit
    R_LG : np array of size (3,3)
        element local2global transformation matrix
    L : float
        element length
    """
    assert 3 == w_G.shape[0]
    lumped_L = np.zeros(12)
    # nodal force
    lumped_L[:3] = w_G*L/2
    lumped_L[6:9] = w_G*L/2
    # transform global load density to local density
    w_L = R_LG.dot(w_G)
    # node 0, 1 local moment
    M_L0 = np.array([0, -w_L[2]*(L**2)/12., w_L[1]*(L**2)/12.])
    M_L1 = -M_L0
    R_GL = R_LG.T
    lumped_L[3:6] = R_GL.dot(M_L0)
    lumped_L[9:12] = R_GL.dot(M_L1)
    return lumped_L

###############################################

class NumpyStiffness(StiffnessBase):
    def __init__(self, nodes, elements, supports, materials, crosssecs, 
        joints=None, verbose=False, output_json=False, unit='meter'):
        assert unit in LENGTH_CONVERSION
        scale = LENGTH_CONVERSION[unit]
        for node in nodes:
            node.point = scale*np.array(node.point)

        super(NumpyStiffness, self).__init__(nodes, elements, 
            supports, materials, crosssecs, joints, verbose, output_json)

        # default displacement tolerance
        self._trans_tol = 1e-3
        self._rot_tol = np.pi/180 * 3
        # gravity settings
        self._include_self_weight_load = True
        # self._full_node_dof = 6
        # self._node_dof = -1
        # internal states
        self._is_init = False
        # element id to dof map
        self._id_map = None
        # node id to dof map
        self._v_id_map = None
        # node id -> element id
        self._node_neighbors = self.get_node_neighbors()
        # elem tag -> elem ind
        self._e_ind_from_tag = self.get_elem_ind_from_elem_tag()

        # (nE x (12x12)) list of element stiffness matrix
        self._element_k_list = []
        # (nE x (12x12)) block-diagnal global2local coordinate transformation matrix
        self._rot_list = []
        # (_ x (6)) dict of lumped element nodal load in global coordinate
        self._element_lumped_nload_list = {}
        # (nE x (12)) list of lumped GRAVITY element nodal load in global coordinate
        self._element_gravity_nload_list = {}
        # (nVx6), np array
        self._nodal_load = None
        # stored computation results
        self._has_stored_deformation = False
        self._stored_existing_ids = []
        self._stored_nD = None # nodal displacement, (nExistNode x 7) np array
        self._stored_eR = None # element reaction, (nExistElement x 13) np array
        self._stored_fR = None # fixities reaction, (nExistFixities x 6) np array
        self._stored_compliance = None
        # output settings
        self._output_json_file_name = ""
        self._output_json_file_path = ""
        # used in non-frame model
        # translational dof id in [0, ..., 2*node_dof-1]
        self._xyz_dof_id = None
        # elemental reaction dof id in [0, ..., 2*node_dof-1]
        self._e_react_dof_id = None

        self._init()

    @property
    def nV(self):
        """number of vertices in the full model
        """
        return len(self._nodes)

    @property
    def nE(self):
        """number of elements in the full model
        """
        return len(self._elements)

    #######################
    # Load input
    def set_self_weight_load(self, include_self_weight=True, gravity_direction=DEFAULT_GRAVITY_DIRECTION):
        self._include_self_weight_load = include_self_weight
        if include_self_weight:
            assert len(gravity_direction) == 3
            gravity_direction = np.array(gravity_direction)
            self_weight_load_density = {}
            for i in range(self.nE):
                crosssec = self.get_element_crosssec(i)
                mat = self.get_element_material(i)
                q_sw = mat.density * crosssec.A
                w_G = gravity_direction * q_sw
                self_weight_load_density[i] = w_G
            self._element_gravity_nload_list = \
                self.compute_element_uniformly_distributed_lumped_load(self_weight_load_density)

    def set_uniformly_distributed_loads(self, unif_dist_loads):
        e_load_density_from_ind = {}
        for e_tag, uload in unif_dist_loads.items():
            for e_ind in self._e_ind_from_tag[e_tag]:
                e_load_density_from_ind[e_ind] = uload.q
        self._element_lumped_nload_list = \
            self.compute_element_uniformly_distributed_lumped_load(e_load_density_from_ind)

    def set_load(self, nodal_forces):
        """[summary]
        
        Parameters
        ----------
        nodal_forces : dict
            {node id : PointLoad}
        """
        self._nodal_load = np.zeros(self.nV*6)
        for v_id, pf in nodal_forces.items():
            self._nodal_load[self._v_id_map[v_id,:]] = np.array(list(pf.force) + list(pf.moment))

    #######################
    # Compute functions
    def set_nodal_displacement_tol(self, transl_tol, rot_tol):
        # in meter, rad
        self._trans_tol = transl_tol
        self._rot_tol = rot_tol

    def compute_dof_permutation(self, exist_element_ids=[]):
        """a row permuatation matrix (multiply left)
        """
        if len(exist_element_ids)==0:
            exist_element_ids = range(self.nE)

        total_dof = 6*self.nV
        # assume all dof does not exist (-1)
        # 0: free, -1: not exist, 1: fixed
        dof_stat = np.ones(total_dof, dtype=int)*-1

        # existing id set
        exist_node_set = set()
        for e in exist_element_ids:
            exist_node_set.update(self._elements[e].end_node_inds)
            # turn the existing node's dof status to free (0)
            dof_stat[self._id_map[e]] = np.zeros(12)

        # assemble the list of fixed dof fixed list, 1x(nFv)
        if len(set(self._supports.keys()) & exist_node_set) == 0:
            if self._verbose : print('Structure floating: no fixed node exists in the partial structure.')
            return False
        n_exist_fixed_nodes = 0
        n_fixed_dof = 0
        for fv_id, support in self._supports.items():
            if fv_id in exist_node_set:
                dof_stat[self._v_id_map[fv_id,:]] = np.array(support.condition, dtype=int)
                # if all existing element bending=False, fix rotational dof
                # to avoid singularity issue
                # TODO partial assembled state?
                # if all([not self._elements[e_id].bending_stiff or e_id not in exist_element_ids for e_id in self._node_neighbors[fv_id]]):
                if all([not self._elements[e_id].bending_stiff for e_id in self._node_neighbors[fv_id]]):
                    # print('pinned dof because of inactive bending')
                    dof_stat[self._v_id_map[fv_id,3:]] = [1,1,1]
                n_exist_fixed_nodes+=1
                assert all([s >= 0 for s in dof_stat[self._v_id_map[fv_id,:]]])
                n_fixed_dof += sum(dof_stat[self._v_id_map[fv_id,:]])

        if n_fixed_dof == 0:
            if self._verbose:
                print('Not stable: At least one node needs to be fixed in the considered substructure!')
            return False
        n_free_dof = 6*len(exist_node_set)-n_fixed_dof

        # compute permutation map
        free_tail = 0
        fix_tail = n_free_dof
        nonexist_tail = 6*len(exist_node_set)
        id_map_RO = np.array(range(total_dof))
        for i in range(total_dof):
            if 0 == dof_stat[i]:
                id_map_RO[free_tail] = i
                free_tail += 1
            elif 1 == dof_stat[i]:
                id_map_RO[fix_tail] = i
                fix_tail += 1
            elif -1 == dof_stat[i]:
                id_map_RO[nonexist_tail] = i
                nonexist_tail += 1
            else:
                raise ValueError

        # a row permuatation matrix (multiply left)
        perm_data = []
        perm_row = []
        perm_col= []
        for i in range(total_dof):
            perm_row.append(i)
            perm_col.append(id_map_RO[i])
            perm_data.append(1)
        Perm = csc_matrix((perm_data, (perm_row, perm_col)), shape=(total_dof, total_dof))
        return Perm, (exist_node_set, n_free_dof, n_fixed_dof, n_exist_fixed_nodes)

    def compute_full_stiffness_matrix(self, exist_element_ids=[]):
        if len(exist_element_ids)==0:
            exist_element_ids = range(self.nE)
        K_full = assemble_global_stiffness_matrix(self._element_k_list, \
            self.nV, self._id_map, exist_element_ids)
        return K_full

    def solve(self, exist_element_ids=[], if_cond_num=False):
        # assume full structure if exist_element_ids is []
        if len(exist_element_ids)==0:
            exist_element_ids = range(self.nE)

        total_dof = 6*self.nV
        
        Perm, (exist_node_set, n_free_dof, n_fixed_dof, n_exist_fixed_nodes) = self.compute_dof_permutation(exist_element_ids)
        PermT = Perm.T

        # permute the full stiffness matrix & carve the needed portion out
        K_full_ = self.compute_full_stiffness_matrix(exist_element_ids)
        K_perm = (Perm.dot(K_full_)).dot(PermT)
        K_mm = K_perm[0:n_free_dof, 0:n_free_dof]
        K_fm = K_perm[n_free_dof:n_free_dof+n_fixed_dof, 0:n_free_dof]

        # TODO make the load vector P sparse
        # assemble load vector
        nodal_load_P_tmp = np.copy(self._nodal_load)
        element_lumped_load = self._create_uniformly_distributed_lumped_load(exist_element_ids)
        nodal_load_P_tmp += element_lumped_load
        if self._include_self_weight_load:
            load_sw = self._create_self_weight_load(exist_element_ids)
            nodal_load_P_tmp += load_sw
        # nodal_load_P_tmp *= -1

        P_m = np.zeros(n_free_dof)
        P_f = np.zeros(n_fixed_dof)
        # reset results
        U_m = np.zeros(n_free_dof)
        self._stored_compliance = 0.0

        # solve
        if norm(nodal_load_P_tmp) > DOUBLE_EPS:
            P_perm = Perm.dot(nodal_load_P_tmp) 
            P_m = P_perm[0:n_free_dof]
            P_f = P_perm[n_free_dof:n_free_dof+n_fixed_dof]

            # TODO Karamba warning when plugging PLA's properties in
            # Material #0 : For isotropic materials 
            # G must be larger than E/3 and smaller than E/2. 
            # If this condition is not fulfilled the material behaves very strange. 
            # It may lead to an unstable structure. May cause errors in exported models.

            # * preconditioning K_mm by K_mm.diagonal()
            # Ax = b
            # (D1 A D2) x_ = (D1) b
            # x = D2 x_
            # take D1 = D2 = [ 1/sqrt(a_{ii}) ]
            if all([abs(a) > DOUBLE_EPS for a in K_mm.diagonal()]):
                # P = diag(1/a_ii) = E_inv.T * E_inv
                # P_inv = diag(a_ii)
                # P_inv = E * E.T
                E_inv_diag = np.array([1/np.sqrt(a) for a in K_mm.diagonal()])
                E_inv = scipy.sparse.diags(E_inv_diag, format='csc')
            else:
                E_inv = scipy.sparse.eye(K_mm.shape[0], format='csc')
            K_mm_precond = (E_inv.dot(K_mm)).dot(E_inv)

            # TODO serialize and check PD and conditioning
            # https://docs.scipy.org/doc/scipy-0.15.1/reference/generated/scipy.sparse.linalg.SuperLU.html#scipy.sparse.linalg.SuperLU
            # * sparse LU decomposition
            # See 'Symmetric or diagonally dominant problems' in: https://portal.nersc.gov/project/sparse/superlu/faq.html
            # tell SuperLU to enter the Symmetric Mode, see example:
            # https://github.com/xiaoyeli/superlu/blob/a3d5233770f0caad4bc4578b46d3b26af99e9c19/EXAMPLE/dlinsol1.c#L60
            try:
                K_LU = SPLA.splu(K_mm_precond, diag_pivot_thresh=0, options=dict(SymmetricMode=True)) 
            except RuntimeError as exception:
                if self._verbose:
                    print ('Sparse solve error: {}'.format(exception))
                
                    eigvals = SPLA.eigsh(K_mm_precond, k=int(K_mm_precond.shape[0]/2.0), which='SM', return_eigenvectors=False) 
                    # eigvals, eigvecs = scipy.linalg.eigh(K_mm_precond.toarray())
                    n_zero_eigs = 0
                    for _, eigval in enumerate(eigvals):
                        if abs(eigval) < DOUBLE_EPS:
                            n_zero_eigs += 1
                            # print(norm(eigvecs[eig_i]))
                    if n_zero_eigs > 0:
                        print('There are {}/{} rigid body modes in the system. This means some parts can move freely without causing deformation.'.format(
                            n_zero_eigs, K_mm.shape[0]))
                    #  Try to use the 'Eigen Modes'-component and activate the display of local coordinate axes: The first eigen-mode will be the rigid body motion.
                    #  If this does not help, check whether you have a pinned support directly attached to a hinge. A hinge introduces an extra node which may cause the problem.
                    #  When analyzing a flat shell structure one has to lock the rotation perpendicular to the plate in at least one node.
                return False

            # * positive definiteness can be checked by
            #   - LDLt, check if D is all positive
            #   - two ends of the spectrum are positive
            # TODO which one is more stable? faster?

            if not (K_LU.U.diagonal()>-DOUBLE_EPS).all():
                # https://gist.github.com/omitakahiro/c49e5168d04438c5b20c921b928f1f5d
                # https://docs.scipy.org/doc/scipy-1.4.1/reference/generated/scipy.sparse.linalg.SuperLU.html#scipy.sparse.linalg.SuperLU
                # Note: the permutation is not checker here
                print('perm: ', (K_LU.perm_r == np.arange(K_mm_precond.shape[0])).all() )
                # when a diagonal entry is smaller than the threshold, the code will still choose an off-diagonal pivot.
                # That is, the row permutation P_r may not be Identity.
                if self._verbose : 
                    print('ERROR: Stiffness Solver fail! Stiffness matrix not Positive Definite, the sub-structure might contain mechanism.')
                return False

            # * checking spectrum positive
            # if False:
            #     try:
            #         # eigen values from both ends of the spectrum
            #         eigvals = SPLA.eigsh(K_mm_precond, k=2, which='BE', tol=1e-2, return_eigenvectors=False) 
            #         if not np.all(eigvals > -DOUBLE_EPS):
            #             if self._verbose:
            #                 print('Stiffness matrix not Positive Definite: largest eigenvalue {}'.format(eigvals))
            #             return False
            #     except SPLA.eigen.arpack.ArpackNoConvergence as errmsg:
            #         # likely a condition number problem
            #         # https://math.stackexchange.com/questions/261295/to-invert-a-matrix-condition-number-should-be-less-than-what
            #         cond_num = np.log(LA.cond(K_mm_precond.toarray()))
            #         if cond_num > -np.log(DOUBLE_EPS):
            #             if self._verbose:
            #                 print('Condition number: {} > 1/eps ({})'.format(np.log(cond_num), -np.log(DOUBLE_EPS)))
            #             return False
            #         if self._verbose:
            #             print('Stiffness matrix condition number explosion: {}'.format(errmsg))
            #         # return False

            # ? Method 1: direct solve
            # U_m_pred = SPLA.spsolve(K_mm_precond, D.dot(P_m))

            # ? Method 2: CG variant
            # https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.gmres.html#scipy.sparse.linalg.gmres
            # U_m, exit_code = SPLA.gmres(K_mm_precond, P_m, tol=1e-7)
            # if exit_code != 0:
            #     if self._verbose:
            #         if exit_code > 0:
            #             print('# iterations: {} - {}'.format(exit_code, 'convergence to tolerance not achieved'))
            #         elif exit_code < 0:
            #             print('exit_code : {} - {}'.format(exit_code, 'illegal input or breakdown'))
            #     return False

            # ? Method 3: SuperLU
            U_m_pred = K_LU.solve(E_inv.dot(P_m))

            # undo the precondition
            U_m = E_inv * U_m_pred

            # gathering results
            self._stored_compliance = 0.5*U_m.dot(P_m)
        else:
            if self._verbose:
                cprint('No external force is applied to the structure.', 'yellow')
            # warnings.warn('No external force is applied to the structure.')

        # reaction
        R = np.zeros(total_dof)
        R[n_free_dof:n_free_dof+n_fixed_dof] = K_fm.dot(U_m) - P_f
        R = PermT.dot(R)
        # displacement
        U = np.zeros(total_dof)
        U[0:n_free_dof] = U_m
        U = PermT.dot(U)

        # raw result conversion        
        self._stored_nD = np.zeros((len(exist_node_set), 1+6), dtype=float)
        for i, v in enumerate(list(exist_node_set)):
            self._stored_nD[i, :] = np.hstack([[v], U[self._v_id_map[v,:]]])

        self._stored_fR = np.zeros((n_exist_fixed_nodes, 1+6), dtype=float)
        cnt = 0
        for fv_id in self._supports:
            if fv_id in exist_node_set:
                self._stored_fR[cnt, :] = np.hstack([[fv_id], R[fv_id*6:(fv_id+1)*6]])
                cnt += 1

        self._stored_eR = np.zeros((len(exist_element_ids), 1+12), dtype=float)
        for i, e in enumerate(exist_element_ids):
            Ue = U[self._id_map[e]]
            eRF = self._rot_list[e].dot(self._element_k_list[e].dot(Ue))
            self._stored_eR[i, :] = np.hstack([[e], eRF])

        self._stored_existing_ids = np.copy(exist_element_ids)
        self._has_stored_deformation = True

        return self.check_stiffness_criteria(self._stored_nD, \
            self._stored_fR, self._stored_eR)

    #######################
    # Public compute functions that are non-state-changing
    def compute_element_uniformly_distributed_lumped_load(self, e_load_density_from_ind):
        """Use fixed-end beam formula to lump uniformly distributed load 
        to equivalent nodal loads
        
        Parameters
        ----------
        element_load_density : dict: int->np array of size (3,)
            load density dict that maps element id to load density vector of size 3
        
        Returns
        -------
        element_lumped_nload_list : dict: e_id : np array of size (3,)
        """
        # refer [MSA McGuire et al.] P111
        # Loads between nodal points
        # TODO avoid allocation here?
        # if isinstance(element_load_density, np.ndarray):
        #     assert element_load_density.shape[1] == 4
        #     tmp_dict = {}
        #     for row in element_load_density:
        #         tmp_dict[int(row[0])] = row[1:]
        #     element_load_density = tmp_dict
        element_lumped_nload_list = {i : np.zeros(12) for i in range(self.nE)}
        for e_ind, w_G in e_load_density_from_ind.items():
            end_u_id, end_v_id = self._elements[e_ind].end_node_inds
            end_u = self._nodes[end_u_id].point 
            end_v = self._nodes[end_v_id].point
            L = norm(end_u-end_v)
            R3 = global2local_transf_matrix(end_u, end_v)
            lumped_L = compute_lumped_uniformly_distributed_load(np.array(w_G), R3, L)
            assert lumped_L.shape[0] == 12
            element_lumped_nload_list[e_ind] = lumped_L
        return element_lumped_nload_list

    def compute_nodal_loads(self, nodal_forces):
        ext_load = np.zeros(6*self.nV)
        for v, nF in nodal_forces.items():
            ext_load[v*6:(v+1)*6] = nF
        return ext_load
    
    def check_stiffness_criteria(self, node_displ, fixities_reaction, element_reaction, ord=None):
        # user can overwrite this function
        nD_trans_max, nD_rot_max, _, _ = self.get_max_nodal_deformation(ord=ord)
        return nD_trans_max < self._trans_tol

    #######################
    # Get functions - Main results
    def has_stored_result(self):
        return self._has_stored_deformation

    def get_solved_results(self):
        assert self.has_stored_result()
        success = self.check_stiffness_criteria(self._stored_nD, self._stored_fR, self._stored_eR)
        return (success, self._stored_nD, self._stored_fR, self._stored_eR)

    def get_compliance(self):
        assert self.has_stored_result()
        return self._stored_compliance

    def get_max_nodal_deformation(self, ord=None):
        assert self.has_stored_result()
        # print(self._stored_nD)
        # print(LA.norm(self._stored_nD[:,1:4], ord=ord, axis=1))
        nD_trans_max = LA.norm(self._stored_nD[:,1:4], ord=ord, axis=1)
        nD_rot_max = LA.norm(self._stored_nD[:,4:7], ord=ord, axis=1)
        return (np.max(nD_trans_max), np.max(nD_rot_max),\
                np.argmax(nD_trans_max), np.argmax(nD_rot_max))

    # element attributes
    def get_element_crosssec(self, elem_id):
        e_tag = self._elements[elem_id].elem_tag
        # assert e_tag in self._crosssecs
        if e_tag in self._crosssecs:
            crosssec = self._crosssecs[e_tag]
        else:
            # TODO: default cross sec, if no [""] key is assigned
            warnings.warn('No cross section assigned for element tag |{}|, using the default tag'.format(e_tag))
            crosssec = self._crosssecs[None]
        return crosssec

    def get_element_material(self, elem_id):
        e_tag = self._elements[elem_id].elem_tag
        # assert e_tag in self._materials
        if e_tag in self._materials:
            mat = self._materials[e_tag]
        else:
            # TODO: default material, if no [""] key is assigned
            warnings.warn('No material assigned for element tag |{}|, using the default tag'.format(e_tag))
            mat = self._materials[None]
        return mat

    # settings getters
    def get_nodal_deformation_tol(self):
        return (self._trans_tol, self._rot_tol)

    def get_frame_stat(self):
        return (self.nV, self.nE)

    #######################
    # Get functions - procedural data
    def get_lumped_nodal_loads(self, existing_ids=[]):
        """packing all applied nodal load and return the load vector
        
        Parameters
        ----------
        existing_ids : list of int, optional
            existing element's indices, by default []
        
        Returns
        -------
        np array of size (nV * 6)

        """
        if len(existing_ids) == 0:
            existing_ids = range(self.nE)
        P = np.copy(self._nodal_load)
        P += self._create_uniformly_distributed_lumped_load(existing_ids)
        if self._include_self_weight_load:
            P += self.get_gravity_nodal_loads(existing_ids)
        return P

    def get_gravity_nodal_loads(self, existing_ids=[]):
        if len(existing_ids)==0:
            existing_ids = range(self.nE)
        return self._create_self_weight_load(existing_ids)

    def get_element_stiffness_matrices(self):
        return self._element_k_list

    def get_element_local2global_rot_matrices(self):
        return self._rot_list

    def get_element2dof_id_map(self):
        return self._id_map

    def get_node2dof_id_map(self):
        return self._v_id_map

    def get_node_neighbors(self):
        """Return all nodes' connected element ids
        
        Returns
        -------
        node_neighbors : dict
            {node_id : set(connected element ids)}
        """
        node_neighbors = defaultdict(set)
        for e in self._elements:
            n1, n2 = e.end_node_inds
            node_neighbors[n1].add(e.elem_ind)
            node_neighbors[n2].add(e.elem_ind)
        return node_neighbors

    def get_elem_ind_from_elem_tag(self):
        e_ind_from_tag = {}
        for e in self._elements:
            if e.elem_tag not in e_ind_from_tag:
                e_ind_from_tag[e.elem_tag] = []
            else:
                e_ind_from_tag[e.elem_tag].append(e.elem_ind)
        return e_ind_from_tag

    ################################
    # Visualizations
    def get_original_shape(self):
        raise NotImplementedError()

    def get_deformed_shape(self, exagg_ratio=1.0, disc=10):
        raise NotImplementedError()

    #######################
    # Output export
    def set_output_json_path(self, file_path, file_name):
        raise NotImplementedError()

    def set_output_json(self, output_json=False):
        raise NotImplementedError()

    ################################
    # private functions - init & precompute
    def _init(self, verbose=False):
        # else:
            # self._node_dof = 6
            # self._xyz_dof_id = np.array(range(2*self._node_dof))
            # self._e_react_dof_id = np.array(range(2*self._node_dof))
        self._precompute_element_stiffness_matrices()
        # node id-to-dof map
        self._v_id_map = np.zeros((self.nV, 6), dtype=int)
        for i in range(self.nV):
            self._v_id_map[i, :] = np.array(range(i*6, (i+1)*6))
        # element id-to-dof map
        self._id_map = np.zeros((self.nE, 12), dtype=int)
        for i in range(self.nE):
            end_u_id, end_v_id = self._elements[i].end_node_inds
            self._id_map[i, :] = np.hstack([self._v_id_map[end_u_id,:], self._v_id_map[end_v_id,:]])
        # loads
        self._nodal_load = np.zeros(self.nV*6)
        if self._include_self_weight_load:
            self.set_self_weight_load(True, DEFAULT_GRAVITY_DIRECTION)
        self._is_init = True
        if verbose:
            print('NpStiffness init done.')

    def _precompute_element_stiffness_matrices(self):
        self._element_k_list = [np.zeros((12,12)) for i in range(self.nE)]
        self._rot_list = [np.zeros((12,12)) for i in range(self.nE)]
        for i in range(self.nE):
            element = self._elements[i]
            e_tag = element.elem_tag
            end_u_id, end_v_id = element.end_node_inds
            end_u = self._nodes[end_u_id].point
            end_v = self._nodes[end_v_id].point
            L = norm(end_u-end_v)

            cr1 = [np.inf for _ in range(3)]
            cr2 = [np.inf for _ in range(3)]
            if e_tag in self._joints:
                cr1 = [np.inf if c is None else c for c in self._joints[e_tag].c_conditions[0]]
                cr2 = [np.inf if c is None else c for c in self._joints[e_tag].c_conditions[1]]
            if not element.bending_stiff:
                cr1 = [0.0 for _ in range(3)]
                cr2 = [0.0 for _ in range(3)]
            R3 = global2local_transf_matrix(end_u, end_v)
            crosssec = self.get_element_crosssec(i)
            mat = self.get_element_material(i)
            K_loc = create_local_stiffness_matrix(L, crosssec.A, crosssec.Jx,\
                crosssec.Iy, crosssec.Iz, mat.E, mat.mu,\
                cr1=cr1, cr2=cr2)
            R_LG = turn_diagblock(R3) # 12x12
            K_G = ((R_LG.T).dot(K_loc)).dot(R_LG)
            self._element_k_list[i] = K_G
            self._rot_list[i] = R_LG

    ################################
    # private functions - solve-time vector/matrix assembly
    def _create_self_weight_load(self, exist_e_ids):
        self_weight_P = np.zeros(self.nV * 6)
        for e in exist_e_ids:
            Qe = self._element_gravity_nload_list[e]
            for j in range(self._id_map.shape[1]):
                self_weight_P[self._id_map[e, j]] += Qe[j]
        return self_weight_P
    
    def _create_uniformly_distributed_lumped_load(self, exist_e_ids):
        element_lump_P = np.zeros(self.nV * 6)
        for e in exist_e_ids:
            if e in self._element_lumped_nload_list:
                Qe = self._element_lumped_nload_list[e]
                for j in range(self._id_map.shape[1]):
                    element_lump_P[self._id_map[e, j]] += Qe[j]
        return element_lump_P