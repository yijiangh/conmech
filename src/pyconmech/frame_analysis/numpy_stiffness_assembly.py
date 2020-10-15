import numpy as np
from numpy.linalg import norm
from scipy.sparse import csc_matrix
from termcolor import cprint

from pyconmech.frame_analysis.io_base import mu2G, G2mu

DOUBLE_EPS = 1e-14
DOUBLE_INF = 1e14

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
    if L < 1e-6:
        raise ValueError('Bar length too small!: Node {} - Node {}'.format(end_vert_u, end_vert_v))

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