import os
import warnings
from pyconmech.frame_analysis.frame_file_io import read_frame_json, read_load_case_json
    
##################################################
# Stiffness Backend Template

class StiffnessBase(object):
    def __init__(self, nodes, elements, supports, materials, crosssecs, 
        joints=None, verbose=False, output_json=False):
        self._nodes = nodes
        self._elements = elements

        # turn lists into dicts
        self._supports = {}
        for support in supports:
            self._supports[support.node_ind] = support

        self._joints = {}
        if joints is not None:
            for joint in joints:
                for e_tag in joint.elem_tags:
                    if e_tag in self._joints:
                        warnings.warn('Multiple joints assigned to the same element tag |{}|!'.format(e_tag))
                    self._joints[e_tag] = joint

        self._materials = {}
        for mat in materials:
            for e_tag in mat.elem_tags:
                if e_tag in self._materials:
                    warnings.warn('Multiple materials assigned to the same element tag |{}|!'.format(e_tag))
                self._materials[e_tag] = mat

        self._crosssecs = {}
        for cs in crosssecs:
            for e_tag in cs.elem_tags:
                if e_tag in self._crosssecs:
                    warnings.warn('Multiple materials assigned to the same element tag |{}|!'.format(e_tag))
                self._crosssecs[e_tag] = cs

        self._verbose = verbose
        self._output_json = output_json
        self._trans_tol = None
        self._rot_tol = None

        if len(self._supports) == 0:
            raise RuntimeError('there needs to be at least one support (fixed) vertex in the model!')

    @classmethod
    def from_json(cls, json_file_path=None, verbose=False):
        """init class from a frame json file
        
        Parameters
        ----------
        json_file_path : [type], optional
            [description], by default None
        verbose : bool, optional
            [description], by default False
        
        Returns
        -------
        [type]
            [description]
        """
        assert os.path.exists(json_file_path), "json file not exists!"
        nodes, elements, supports, joints, materials, crosssecs, _, _ = \
            read_frame_json(json_file_path, verbose=verbose)
        return cls(nodes, elements, supports, materials, crosssecs, 
            joints=joints, verbose=verbose)

    #######################
    # Load input

    def set_self_weight_load(self, include_self_weight=True, gravity_direction=[0,0,-1]):
        raise NotImplementedError()

    def set_load(self, nodal_forces):
        raise NotImplementedError()

    def set_uniformly_distributed_loads(self, element_load_density):
        raise NotImplementedError()

    def _parse_load_case_from_json(self, file_path):
        # TODO cpp engine only returns nodal load now
        point_load, uniform_element_load, include_self_weight = \
            read_load_case_json(file_path)
        return point_load, uniform_element_load, include_self_weight

    #######################
    # Output write

    def set_nodal_displacement_tol(self, transl_tol, rot_tol):
        raise NotImplementedError()

    def solve(self, exist_element_ids=[], if_cond_num=False):
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
        raise NotImplementedError()

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

