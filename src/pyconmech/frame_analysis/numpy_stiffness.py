"""
# Variables:
displacements:          {u}
strain-displacement eq: {\epsilon} = [\partial]{u}

# Induced field:
constituitive matrix:    elastic constants [E]
stress-strain eq:        stresses {\sigma} = [E]{\epsilon}

# Boundary condition:
surface load:         traction {\Phi} on boundary
internal load:        body force F_x, F_y, F_z

# Compatibility:
displacement continuous across the body

equilibrium eq:       [\partial]^T {\sigma} + {F} = 0

# principle of virtual work:
For any *quasi-static* and admissible virtual displacement {\delta u} from
an equilibrium conf, the increment of strain energy stored is equal to
the increment of work done by body (internal) forces {F} in volume V and
traction force {\Phi} on boundary S.

    \int {\delta \epsilon}^T {\sigma} dV = 
        \int {\delta u}^T F dV + \int {\delta u}^T {\Phi} dS
"""

import numpy as np
from numpy.linalg import norm, solve
from scipy.sparse import csr_matrix
from numpy.polynomial.polynomial import polyvander, polyval, Polynomial
from .stiffness_base import StiffnessBase

DOUBLE_EPS = 1e-14

# internal calculation unit
LENGTH_UNIT = 'meter'
ROT_ANGLE_UNIT = 'rad'
FORCE_UNIT = 'kN'

# derived units
MOMENT_UNIT = FORCE_UNIT + "*" + ROT_ANGLE_UNIT
AREA_UNIT = LENGTH_UNIT + "^2"
AREA_INERTIA_UNIT = LENGTH_UNIT + "^4"
PRESSURE_UNIT = FORCE_UNIT + "/" + LENGTH_UNIT + "^2"; # "kN/m^2" modulus
DENSITY_UNIT  = FORCE_UNIT + "/" + LENGTH_UNIT + "^3"; # "kN/m^3"


def axial_stiffness_matrix(L, A, E):
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

    Return
    ------
    K_ax_x : 2x2 numpy array
    """
    K = np.ones([2,2])
    K[0,1] = -1.
    K[1,0] = -1.
    K *= (A*E/L)
    return K


def torsional_stiffness_matrix(L, J, G):
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

    Return
    ------
    K_tor_x : 2x2 numpy array
    """
    return axial_stiffness_matrix(L, J, G)


def bending_stiffness_matrix(L, E, Iz, axis=2):
    """stiffness matrix of uniform beam element without the 
    transverse shear deformation

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

    Return
    ------
    K_bend_z : 4x4 numpy array
    """
    K = np.zeros([4,4])
    sign = 1. if axis==2 else -1.
    K[0,:] = np.array([12.,      sign*6*L, -12.,      sign*6*L])
    K[1,:] = np.array([sign*6*L, 4*(L**2), sign*-6*L, 2*(L**2)])
    K[2,:] = -K[0,:]
    K[3,:] = np.array([sign*6*L, 2*(L**2), -sign*6*L, 4*(L**2)])
    K *= (E*Iz/L**3)
    return K

def mu2G(mu, E):
    return E/(2*(1+mu))

def G2mu(G, E):
    return E/(2*G)-1

def create_local_stiffness_matrix(L, A, Jx, Iy, Iz, E, mu):
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
    
    Returns
    -------
    K : 12x12 numpy array
    """
    G = mu2G(mu, E)
    # Fx1, Fx2 : u1, u2
    # 0, 6 : 0, 6
    axial_x_k = axial_stiffness_matrix(L, A, E)
    # Mx1, Mx2 : \theta_x1, \theta_x2
    # 3, 9 : 3, 9
    tor_x_k = torsional_stiffness_matrix(L, Jx, G)

    # Fy1, Mz1, Fy2, Mz2 : v1, \theta_z1, v2, \theta_z2
    # 1, 5, 7, 11 : 1, 5, 7, 11
    bend_z_k = bending_stiffness_matrix(L, E, Iz, axis=2)
    # Fz1, My1, Fz2, My2 : v1, \theta_z1, v2, \theta_z2
    # 2, 4, 8, 10 : 2, 4, 8, 10
    bend_y_k = bending_stiffness_matrix(L, E, Iz, axis=1)

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
    Ksp : scipy.sparse.csr_matrix
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
    Ksp = csr_matrix((data, (row, col)), shape=(total_dof, total_dof))
    return Ksp

# LINEAR_INTERP_POLY = (np.array([1, -1/L]),
#                       np.array([0, 1/L]))
# CUBIC_INTERP_POLY = (
#     np.array([1, 0, -3/L**2, 2/L**3]),
#     np.array([0, 1, -2/L, 1/L**2]),
#     np.array([0, 0, 3/L**2, 2/L**3]),
#     np.array([0, 0, -1/L, 1/L**2]))

def get_element_shape_fn(end_u, end_v, d_u, d_v, exagg=1.0):
    """cubic polynomial interpolation given its boundary condition:
        d = [u1, du1/dx, u2, du2/dx]
    
    Parameters
    ----------
    d : [type]
        [description]
    L : float
        element length

    Return
    ------
    poly_eval_fn : function handle
    """
    L = norm(np.array(end_u) - np.array(end_v))

    R3 = global2local_transf_matrix(end_u, end_v)
    R = np.zeros((12,12))
    for i in range(4):
        R[i*3:(i+1)*3, i*3:(i+1)*3] = R3

    # compute end deflection in local ordinates
    D_local = np.zeros((12, 1))
    D_global = np.hstack([d_u, d_v])
    D_local = exagg * R.dot(D_global)

    u_poly = []
    for i in range(3):
        if i==0:
            # linear monomials 1+x
            A = polyvander([0, L], 1)
            coeff = solve(A, [D_local[i], D_local[i+6]])
        else:
            # cubic monomials
            A = np.zeros((4, 4))
            A[[0,2], :] = polyvander([0, L], 3)
            c = np.hstack([Polynomial([1]*4).deriv().coef, [0]])
            A[[1,3], :] = polyvander([0, L], 3) * c

            # TODO careful about the sign
            coeff = solve(A, [D_local[i], D_local[i+6]])
        u_poly.append(coeff)
    
    def shape_fn(t):
        assert t>=0 and t<=1.0
        d = np.zeros(3)
        # TODO: transf back to global
        for i in range(3):
            # d[i] = 
        # d = np.array([np.polynomial.polynomial.polyval(t, u_poly[i]) for i in range(3)])
        return end_u + d

    return shape_fn

class NumpyStiffness(StiffnessBase):
    def __init__(self, vertices, elements, fixities, material_dicts, 
        verbose=False, model_type='frame', output_json=False):
        super(NumpyStiffness, self).__init__(vertices, elements, fixities, material_dicts, \
            verbose, model_type, output_json)
        # default displacement tolerance
        self._trans_tol = 1e-3
        self._rot_tol = np.pi/180 * 3

    #######################
    # Load input

    def set_self_weight_load(self, gravity_direction, include_self_weight=True):
        raise NotImplementedError()

    def set_load(self, nodal_forces):
        """input: #nL x 7 numpy matrix
        """
        raise NotImplementedError()

    def set_uniformly_distributed_loads(self, element_load_density):
        raise NotImplementedError()

    def _parse_load_case_from_json(self, file_path):
        raise NotImplementedError()

    #######################
    # Output write

    def set_nodal_displacement_tol(self, transl_tol, rot_tol):
        raise NotImplementedError()

    def solve(self, exist_element_ids=[], if_cond_num=True):
        raise NotImplementedError()

    #######################
    # Output export

    def set_output_json_path(self, file_path, file_name):
    # .def("set_output_json_path", &conmech::stiffness_checker::Stiffness::setOutputJsonPath,
    #   py::arg("file_path"), py::arg("file_name"))
        raise NotImplementedError()

    def set_output_json(self, output_json=False):
    # .def("set_output_json", &conmech::stiffness_checker::Stiffness::setOutputJson,
    #   py::arg("output_json") = false)
        raise NotImplementedError()

    #######################
    # Get functions

    def get_frame_stat(self):
        return 

    def get_lumped_nodal_loads(self, existing_ids=[]):
        raise NotImplementedError()

    def get_gravity_nodal_loads(self, existing_ids=[]):
        raise NotImplementedError()

    def get_element_stiffness_matrices(self):
        raise NotImplementedError()

    def get_element_local2global_rot_matrices(self):
        raise NotImplementedError()

    def get_element2dof_id_map(self):
        raise NotImplementedError()

    def get_node2dof_id_map(self):
        raise NotImplementedError()

    def has_stored_result(self):
        raise NotImplementedError()

    def get_solved_results(self):
        raise NotImplementedError()

    def get_max_nodal_deformation(self):
        raise NotImplementedError()

    def get_nodal_deformation_tol(self):
        raise NotImplementedError()

    def get_original_shape(self):
        raise NotImplementedError()

    def get_deformed_shape(self, exagg_ratio=1.0, disc=10):
        raise NotImplementedError()

    ################################
    # private functions
    def _create_complete_global_stiffness_matrix(self):
        pass