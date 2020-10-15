import warnings
import numpy as np
import scipy
import scipy.sparse.linalg as SPLA
from scipy.sparse import csc_matrix
import numpy.linalg as LA
from collections import defaultdict
from numpy.linalg import norm
from termcolor import cprint

from pyconmech.frame_analysis.stiffness_base import StiffnessBase
from pyconmech.frame_analysis.numpy_stiffness_assembly import DOUBLE_EPS, DOUBLE_INF, \
    create_local_stiffness_matrix, global2local_transf_matrix, turn_diagblock, assemble_global_stiffness_matrix, \
    compute_lumped_uniformly_distributed_load

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
        self._include_self_weight_load = False
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
        if unif_dist_loads is not None:
            e_load_density_from_ind = {}
            for e_tag, uload in unif_dist_loads.items():
                for e_ind in self._e_ind_from_tag[e_tag]:
                    e_load_density_from_ind[e_ind] = uload.q
            self._element_lumped_nload_list = \
                self.compute_element_uniformly_distributed_lumped_load(e_load_density_from_ind)
        else:
            # reset
            self._element_lumped_nload_list = {}

    def set_load(self, nodal_forces):
        """[summary]
        
        Parameters
        ----------
        nodal_forces : dict
            {node id : PointLoad}
        """
        if nodal_forces is not None:
            for v_id, pf in nodal_forces.items():
                self._nodal_load[self._v_id_map[v_id,:]] = np.array(list(pf.force) + list(pf.moment))
        else:
            # reset
            self._nodal_load = np.zeros(self.nV*6)

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
        if len(self._element_lumped_nload_list) > 0:
            # if element unif load is specified
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

            # ? Method 1: direct solve
            # U_m_pred = SPLA.spsolve(K_mm_precond, E_inv.dot(P_m))

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
        element_lumped_nload_list = {i : np.zeros(12) for i in range(self.nE)}
        # element_lumped_nload_list = {}
        for e_ind, w_G in e_load_density_from_ind.items():
            end_u_id, end_v_id = self._elements[e_ind].end_node_inds
            end_u = self._nodes[end_u_id].point 
            end_v = self._nodes[end_v_id].point
            L = norm(end_u-end_v)
            try:
                R3 = global2local_transf_matrix(end_u, end_v)
            except ValueError as err:
                raise ValueError('E#{} - {}'.format(e_ind, err))
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
        if len(self._element_lumped_nload_list) > 0:
            # if element unif load is specified
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
            try:
                R3 = global2local_transf_matrix(end_u, end_v)
            except ValueError as err:
                raise ValueError('E#{} - {}'.format(i, err))
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
        assert len(self._element_gravity_nload_list) > 0, 'Call `set_self_weight_load` first to init the load!'
        self_weight_P = np.zeros(self.nV * 6)
        for e in exist_e_ids:
            Qe = self._element_gravity_nload_list[e]
            for dof_j in range(self._id_map.shape[1]):
                self_weight_P[self._id_map[e, dof_j]] += Qe[dof_j]
        return self_weight_P
    
    def _create_uniformly_distributed_lumped_load(self, exist_e_ids):
        assert len(self._element_lumped_nload_list) > 0, 'Call `set_load` first to init the load!'
        element_lump_P = np.zeros(self.nV * 6)
        for e in exist_e_ids:
            Qe = self._element_lumped_nload_list[e]
            for dof_j in range(self._id_map.shape[1]):
                element_lump_P[self._id_map[e, dof_j]] += Qe[dof_j]
        return element_lump_P