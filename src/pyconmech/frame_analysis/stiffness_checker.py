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

    def set_loads(self, point_loads={}, include_self_weight=False, uniform_distributed_load={}):
        self._sc_ins.set_self_weight_load(include_self_weight)
        if point_loads:
            pt_loads = []
            for vid, vload in point_loads.items():
                pt_loads.append([vid] + vload)
            self._sc_ins.set_load(np.array(pt_loads))
        if uniform_distributed_load:
            raise NotImplementedError

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
    
    def get_solved_results(self):
        """Fetch back solved results from last time.

        The elemental local coordinate is placed at the end point 0, with local axis pointing along
        the `end point 0` -> `end point L` direction. See doc (TODO: link) for more info.
        
        Returns
        -------
        success : bool
            pass criteria or not
        nD : numpy array
            nodal displacements in the global coordinate
            [[node_id, dx, dy, dz, rxx, ryy, rzz], ...]
        fR : numpy array
            fixities reaction force and moment in the global coordinate
            [[node_id, Fx, Fy, Fz, Mxx, Myy, Mzz], ...]
        eR : numpy array
            elemental reaction force and moment in the local coordinate
            [[element_id, F_0_lx, F_0_ly, F_0_lz, M_0_lxx, M_0_lyy, M_0_lzz,
                          F_L_lx, F_L_ly, F_L_lz, M_L_lxx, M_L_lyy, M_L_lzz,
            ], ...]
            F_0_lx means internal reaction force at the end point 0, in the direction of local x axis
            M_L_lyy means internal reaction moment at the end point L, around the local yy axis            
        """
        success, nD, fR, eR = self._sc_ins.get_solved_results()
        return success, nD, fR, eR


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


