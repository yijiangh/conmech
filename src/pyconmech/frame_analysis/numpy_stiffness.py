import numpy as np
from .stiffness_base import StiffnessBase

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
    K_ax_x : 2x2 numpy matrix
    """
    K = np.ones([2,2])
    K[0,1] = -1.
    K[1,0] = -1.
    return K * (A*E/L)


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
    K_tor_x : 2x2 numpy matrix
    """
    return axial_stiffness_matrix(L, J, G)


def bending_stiffness_matrix(L, E, Iz):
    """uniform beam element without transverse shear deformation

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

    Return
    ------
    K_bend_z : 4x4 numpy matrix
    """
    K = np.zeros([4,4])
    K[0,:] = np.array([12., 6*L, -12., 6*L])
    K[1,:] = np.array([6*L, 4*(L**2), -6*L, 2*(L**2)])
    K[2,:] = -K[0,:]
    K[3,:] = np.array([6*L, 2*(L**2), -6*L, 4*(L**2)])
    return K*(E*Iz/L**3)


def create_local_stiffness_matrix(L, A, Jx, Iy, Iz, E, G, mu, dim=3):
    # if dim != 3:
    raise NotImplementedError()
    # assembly components

def global2local_transf_matrix():
    pass

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