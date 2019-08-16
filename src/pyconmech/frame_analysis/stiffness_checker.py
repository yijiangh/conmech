"""
Python interface class for the stiffness checker backend.
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from tempfile import TemporaryDirectory

import os
from collections import defaultdict
import numpy as np

from _pystiffness_checker import _stiffness_checker
from pyconmech.frame_analysis.frame_file_io import read_frame_json, write_frame_json

class stiffness_checker(object):
    """stiffness checking instance for 3D frame deformation analysis

    Calculating elastic deformation given a 3D frame shape specified in a json file.
    See tests/test_data for examples on the frame shape json file format.
    
    """
    def __init__(self, json_file_path, verbose=False):
        """Init fn for stiffness_checker

        By default, self-weight load is added.
        Disable by <your stiffness checker class>.set_self_weight_load(False)
        
        Parameters
        ----------
        json_file_path : str
            absolute path to the frame shape's json file
        verbose : bool, optional
            verbose screen outputs turned on/off, by default False
        """
        assert os.path.exists(json_file_path)
        self._sc_ins = _stiffness_checker(json_file_path=json_file_path, verbose=verbose)
        self._sc_ins.set_self_weight_load(True)
        
        node_points, element_vids, fix_node_ids, fix_specs, model_type, material_dict, model_name = \
        read_frame_json(json_file_path, verbose=verbose)
        if model_type == 'truss':
            raise NotImplementedError('truss model is not supported now...')

        # TODO: wrap them into a frame class, and provide __eq__ and __ne__
        self._model_name = model_name
        self._node_points = node_points
        self._elements = element_vids
        self._fix_node_ids = fix_node_ids
        self._fix_specs = fix_specs
        self._model_type = model_type
        self._material_dict = material_dict
        
        self.set_nodal_displacement_tol()

    @classmethod
    def from_json(cls, json_file_path, verbose=False):
        return cls(json_file_path, verbose=verbose)

    @classmethod
    def from_frame_data(cls, nodes, elements, fixed_node_ids, material_dict, 
        fixity_specs={}, unit='meter', model_type='frame', model_name=None, verbose=False):

        # here = os.path.dirname(os.path.abspath(__file__))
        # tmp_path = os.path.join(here, 'pyconmech_frame_temp.json')
        with TemporaryDirectory() as temp_dir:
            tmp_path = os.path.join(temp_dir, 'pyconmech_frame_temp.json')
            write_frame_json(tmp_path, nodes, elements, fixed_node_ids, material_dict, 
            fixity_specs=fixity_specs, model_type=model_type, model_name=model_name)
            instance = cls.from_json(tmp_path, verbose)
        return instance

    # ==========================================================================
    # properties
    # ==========================================================================

    @property
    def model_name(self):
        return self._model_name

    @property
    def node_points(self):
        return self._node_points
    
    @property
    def elements(self):
        return self._elements
    
    @property
    def fix_node_ids(self):
        return self._fix_node_ids

    @property
    def fix_element_ids(self):
        fix_e_ids = []
        for e_id, e in enumerate(self.elements):
            if any([v_id in self.fix_node_ids for v_id in e]):
                fix_e_ids.append(e_id)
        return fix_e_ids

    @property
    def fix_specs(self):
        return self._fix_specs

    @property
    def model_type(self):
        return self._model_type

    @property
    def material(self):
        return self._material_dict

    # ==========================================================================
    # solve functions
    # ==========================================================================

    def solve(self, exist_element_ids=[], if_cond_num=False, eid_sanity_check=False):
        """compute elastic deformation using the linear stiffness equation, return
        True or False indicating if the set criteria is passed nor not.
        
        Parameters
        ----------
        exist_element_ids : list, optional
            list of existing element ids, by default []
            By default (exist_element_ids = []), the stiffness checker will 
            compute for the full structure.
        if_cond_num : bool, optional
            check condition number of the stiffness matrix, by default False
            This will increase the fidelity of the computation but induce some overhead.
            Condition number checking is currently not implemented.
        eid_sanity_check : bool, optional
            check if the given element ids are within range. This is simply a sanity check
            but might induce some extra overhead if the solve function is called many times.
            by default False.
        
        Returns
        -------
        bool
            success or not
        """
        if eid_sanity_check:
            e_id_range = list(range(len(self.elements)))
            for e_id in exist_element_ids:
                assert e_id in e_id_range, 'element id not within range!'
        return self._sc_ins.solve(exist_element_ids, if_cond_num)

    # ==========================================================================
    # load settings
    # ==========================================================================

    def set_loads(self, point_loads=None, include_self_weight=False, uniform_distributed_load={}):
        """set load case for the stiffness checker.
        
        Parameters
        ----------
        point_loads : dict, optional
            {node_id : [Fx, Fy, Fz, Mxx, Myy, Mzz]}, in global coordinate,
            by default {}
        include_self_weight : bool, optional
            include gravity load or not, by default False
        uniform_distributed_load : dict, optional
            elemental uniformly distributed load, by default {}
        
        Raises
        ------
        NotImplementedError
            uniform loads input not supported for now
        """
        self._sc_ins.set_self_weight_load(include_self_weight)
        if point_loads:
            pt_loads = []
            for vid, vload in point_loads.items():
                pt_loads.append([vid] + vload)
            self._sc_ins.set_load(np.array(pt_loads))
        if uniform_distributed_load:
            raise NotImplementedError('uniform loads input not supported for now')

    def set_self_weight_load(self, include_self_weight):
        """Turn on/off self-weight load.
        
        Parameters
        ----------
        include_self_weight : bool
        """
        self._sc_ins.set_self_weight_load(include_self_weight)

    # ==========================================================================
    # check criteria settings
    # ==========================================================================

    def set_nodal_displacement_tol(self, trans_tol=1e-3, rot_tol=0.1745):
        """Set nodal displacement tolerance for stiffness checking criteria.

        If the maximal nodal displacement after the stiffness solve exceeds the tolerance,
        then stiffness_checker.solve() (or with input [ids]) will return False.
        
        Parameters
        ----------
        trans_tol : float, optional
            maximal allowable nodal translational movement in global coordinate, unit in meter, by default 1e-3
        rot_tol : float, optional
            maximal allowable nodal rotational movement in global coordinate, unit in meter, by default 0.1745 
            (about 10 degrees)
        """
        self._sc_ins.set_nodal_displacement_tol(trans_tol, rot_tol)

    def get_nodal_deformation_tol(self):
        """Get nodal displacement tolerance for stiffness checking criteria.
        
        Returns
        -------
        trans_tol : float
            maximal nodal translational deformation, entry-wise
        rot_tol : float
            maximal nodal rotational deformation, entry-wise
        """
        trans_tol, rot_tol = self._sc_ins.get_nodal_deformation_tol()
        return trans_tol, rot_tol

    # TODO: pass in user-specified criteria checking function handle

    # ==========================================================================
    # stiffness data query
    # ==========================================================================

    def has_stored_result(self):
        """Check if the stiffness checker has solved results stored,
        i.e. solve function has been called at least once

        Returns
        -------
        bool
        """
        return self._sc_ins.has_stored_result()
    
    def get_solved_results(self, existing_ids=[]):
        """Fetch back solved results from last time.
        
        Returns
        -------
        success : bool
            pass criteria or not
        nD : dict
            nodal displacements in the global coordinate
            {node_id: np.array([dx, dy, dz, rxx, ryy, rzz]), ...]
        fR : dict
            fixities reaction force and moment in the global coordinate
            {node_id: np.array([Fx, Fy, Fz, Mxx, Myy, Mzz], ...}
        eR : numpy array
            elemental reaction force and moment in the local coordinate
            {element_id: {0 : np.array([F_0_lx, F_0_ly, F_0_lz, M_0_lxx, M_0_lyy, M_0_lzz]),
                          1 : np.array([F_L_lx, F_L_ly, F_L_lz, M_L_lxx, M_L_lyy, M_L_lzz])
            ], ...]
            F_0_lx means internal reaction force at the end point 0, in the direction of local x axis
            M_L_lyy means internal reaction moment at the end point L, around the local yy axis            

            The elemental local coordinate is placed at the end point 0, with local axis pointing along
            the `end point 0` -> `end point L` direction. See doc (TODO: link) for more info.
        """
        if not existing_ids:
            existing_e_ids = list(range(len(self.elements)))
        existing_n_ids = self.get_element_connected_node_ids(existing_ids=existing_ids)
        success, nD_np, fR_np, eR_np = self._sc_ins.get_solved_results()
        nD = {}
        for row in nD_np:
            if int(row[0]) in existing_n_ids:
                nD[int(row[0])] = np.array(row[1:7])
        fR = {}
        for row in fR_np:
            if int(row[0]) in existing_n_ids and int(row[0]) in self.fix_node_ids:
                fR[int(row[0])] = np.array(row[1:7])
        eR = {}
        for row in eR_np:
            if int(row[0]) in existing_e_ids:
                e_r = {}
                e_r[0] = np.array(row[1:7])
                e_r[1] = np.array(row[7:13])
                eR[int(row[0])] = e_r
        return success, nD, fR, eR

    def get_max_nodal_deformation(self):
        """Get max nodal deformation info
        
        Returns
        -------
        max_trans : float
            maximal translational deformation, componentwise max, unit meter
        max_rot : float
            maximal rotational deformation, componentwise max, unit rad
        max_trans_vid : int
            node id for maximal trans deformation
        max_rot_vid : int
            node id for maximal rotational deformation
        """
        max_trans, max_rot, max_trans_vid, max_rot_vid = self._sc_ins.get_max_nodal_deformation()
        return max_trans, max_rot, max_trans_vid, max_rot_vid

    def get_compliance(self):
        """Get compliance of the last solved deformation

        compliance = nodal_loads.dot(nodal_deformation)  
        TODO: check if we need moment and rotational deformation here?
        
        Returns
        -------
        float
        """
        return self._sc_ins.get_compliance()

    def get_self_weight_loads(self, existing_ids=[], dof_flattened=False):
        """Return lumped gravity loads

        TODO: put a link to doc here
        
        Parameters
        ----------
        existing_ids : list, optional
            existing element's ids in the partial structure, by default [], which means full structure
        dof_flattened : bool, optional
            if True, return a flattened dof x 1 vector, 
            otherwise return a dict{node_id: np.array([1:dof])}, by default False
        
        Returns
        -------
        dict (dof_flattened = False)
            {node_id : [Fx, Fy, Fz, Mxx, Myy, Mzz]} in global coordinate
        or np.array (dof_flattened = True)
            [Fx, Fy, Fz, Mxx, Myy, Mzz, ...]
        """
        assert self.model_type == 'frame', 'this function assumes 6 dof each node for now!'
        nL_flat = self._sc_ins.get_gravity_nodal_loads(existing_ids)
        if not dof_flattened:
            nL = {}
            for nid in range(len(self.node_points)):
                nL[nid] = nL_flat[nid*6 : (nid+1)*6]
            return nL
        else:
            return nL_flat


    def get_nodal_loads(self, existing_ids=[], dof_flattened=False):
        """Return nodal loads
        
        Parameters
        ----------
        existing_ids : list, optional
            existing element's ids in the partial structure, by default [], which means full structure
        dof_flattened : bool, optional
            if True, return a flattened dof x 1 vector, 
            otherwise return a dict{node_id: np.array([1:dof])}, by default False
        
        Returns
        -------
        dict (dof_flattened = False)
            {node_id : [Fx, Fy, Fz, Mxx, Myy, Mzz]} in global coordinate
        or np.array (dof_flattened = True)
            [Fx, Fy, Fz, Mxx, Myy, Mzz, ...]
        """
        assert self.model_type == 'frame', 'this function assumes 6 dof each node for now!'
        nL_flat = self._sc_ins.get_lumped_nodal_loads(existing_ids)
        if not dof_flattened:
            nL = {}
            for nid in range(len(self.node_points)):
                nL[nid] = nL_flat[nid*6 : (nid+1)*6]
            return nL
        else:
            return nL_flat

    def get_element_stiffness_matrices(self, in_local_coordinate=False):
        """Get all elemental stiffness matrices, each of which is 12 x 12
        in global coordinate (default): R_{LG}^T * K_{eL} * R_{LG}
        in local coordinate (default): K_{eL}
        
        Returns
        -------
        dict
            {element_id : np.array}
        """
        eMs = {}
        eMs_raw = self._sc_ins.get_element_stiffness_matrices()
        if in_local_coordinate:
            R_lgs = self.get_element_local2global_rot_matrices()

        for eid, eMg in enumerate(eMs_raw):
            if not in_local_coordinate:
                eMs[eid] = eMg
            else:
                R_lg = R_lgs[eid]
                eMs = np.transpose(R_lg).dot(eMg)
                eMs = eMs.dot(R_lg)
                eMs[eid] = eMg
        return eMs

    def get_element_local2global_rot_matrices(self):
        """Get all elemental local to global transformation matrices, each of which is 12 x 12

        The array is in the following shape:
        | R       |
        |   R     |
        |     R   |
        |       R |
        where R is the 3x3 coordinate transformation matrix, which tranforms 
        the global xyz axis to the element local axis.

        The coordinate transformation matrix can be used to:
         - transform frame element end forces from the element (local) coordinate system
           to the structure (global) coordinate system
         - transfrom end displacements from the structural (global) coordinate system 
           to the element (local) coordinate system,
         - transform the frame element stiffness and mass matrices
           from element (local) coordinates to structral (global) coordinates.
        Symbolically, the return matrix R = {local}_R_{global}

        Returns
        -------
        dict
            {node_id : 12 x 12 np.array}
        """
        return {nid : mat for nid, mat in enumerate(self._sc_ins.get_element_local2global_rot_matrices())}

    def get_element2dof_id_map(self):
        """Get element_id-to-dof_id maps
        
        Returns
        -------
        dict
            {e_id : {0 : [dof_ids for end point 0]}, {1 : [dof_ids for end point 1]}}
        """
        assert self.model_type == 'frame', 'this function assumes 6 dof each node for now!'
        return {int(e_id) : {0 : id_map[0:6], 1 : id_map[6:12]} 
                for e_id, id_map in enumerate(self._sc_ins.get_element2dof_id_map())}

    def get_node2dof_id_map(self):
        """Get node_id-to-dof_id maps
        
        Returns
        -------
        dict
            {node_id : [dof ids]}
        """
        return {int(n_id) : id_map for n_id, id_map in enumerate(self._sc_ins.get_node2dof_id_map())}

    # ==========================================================================
    # output settings
    # ==========================================================================

    def set_output_json(self, do_output=True, output_dir=None, file_name=None):
        self._sc_ins.set_output_json(do_output)
        if do_output:
            if not output_dir:
                output_dir = os.path.dirname(os.path.abspath(__file__))
                print('NOTE: no output path is specified, the result file is saved at:\n {}'.format(output_dir))
            else:
                if output_dir[-1] != os.sep and output_dir[-2:] != os.sep:
                    output_dir += os.sep # path sep required by the cpp extmodule
            if not file_name:
                file_name = self.model_name + '_result' + '.json'
            self._sc_ins.set_output_json_path(output_dir, file_name)
    
    # ==========================================================================
    # frame data query (TODO: moved to frame class)
    # ==========================================================================

    def get_node_neighbors(self):
        """Return all nodes' connected element ids
        
        Returns
        -------
        node_neighbors : dict
            {node_id : set(connected_element_ids)}
        """
        node_neighbors = defaultdict(set)
        for e in self.elements:
            n1, n2 = e
            node_neighbors[n1].add(e)
            node_neighbors[n2].add(e)
        return node_neighbors
    
    def get_element_neighbors(self):
        """Return all elements' connected element ids
        
        Returns
        -------
        element_neighbors : dict
            {element_id : set(connected element_ids)}
        """
        node_neighbors = self.get_node_neighbors()
        element_neighbors = defaultdict(set)
        for e in self.elements:
            n1, n2 = e
            element_neighbors[e].update(node_neighbors[n1])
            element_neighbors[e].update(node_neighbors[n2])
            element_neighbors[e].remove(e)
        return element_neighbors

    def get_element_connected_node_ids(self, existing_ids=[], fix_node_only=False):
        if not existing_ids:
            existing_ids = list(range(len(self.elements)))
        connected_node_ids = set()
        for e_id in existing_ids:
            if not fix_node_only:
                connected_node_ids.update([v_id for v_id in self.elements[e_id]])
            else:
                connected_node_ids.update([v_id for v_id in self.elements[e_id] if v_id in self.fix_node_ids])
        return list(connected_node_ids)

    # ==========================================================================
    # beam shape query for deformation visualization
    # ==========================================================================

    def get_original_shape(self, disc=1, draw_full_shape=True):
        pass

    def get_deformed_shape(self, disc=10, exagg_ratio=1.0):
        pass