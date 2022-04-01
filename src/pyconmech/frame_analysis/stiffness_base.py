import os
import warnings
from pyconmech.frame_analysis.io_base import Model, LoadCase
    
##################################################
# Stiffness Backend Template

class StiffnessBase(object):
    def __init__(self, model, verbose=False, output_json=False):
        self._nodes = model.nodes
        self._elements = model.elements
        self._supports = model.supports # dict indexed by node_ind
        self._joints = model.joints
        self._joint_id_from_etag = model.joint_id_from_etag
        self._materials = model.materials
        self._material_id_from_etag = model.material_id_from_etag
        self._crosssecs = model.crosssecs
        self._crosssec_id_from_etag = model.crosssec_id_from_etag

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
        return cls(Model.from_json(json_file_path, verbose=verbose), verbose=verbose)

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
        raise NotImplementedError()
        # lc = LoadCase.from_json(file_path)
        # point_load, uniform_element_load, include_self_weight = \
        #     read_load_case_json(file_path)
        # return point_load, uniform_element_load, include_self_weight

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
        assert e_tag in self._crosssec_id_from_etag
        if e_tag in self._crosssec_id_from_etag:
            cs_id = self._crosssec_id_from_etag[e_tag]
        # else:
        #     # warnings.warn('No cross section assigned for element tag |{}|, using the default tag'.format(e_tag))
        #     cs_id = self._crosssec_id_from_etag[None] if None in self._crosssec_id_from_etag else self._crosssec_id_from_etag['']
        return self._crosssecs[cs_id]

    def get_element_material(self, elem_id):
        e_tag = self._elements[elem_id].elem_tag
        assert e_tag in self._material_id_from_etag
        if e_tag in self._material_id_from_etag:
            m_id = self._material_id_from_etag[e_tag]
        # else:
            # TODO: default material, if no key is found
            # warnings.warn('No material assigned for element tag |{}|, using the default tag'.format(e_tag))
            # m_id = self._material_id_from_etag[None] if None in self._material_id_from_etag else self._material_id_from_etag['']
        return self._materials[m_id]

    def get_element_joint(self, elem_id):
        e_tag = self._elements[elem_id].elem_tag
        if e_tag in self._joint_id_from_etag:
            jt_id = self._joint_id_from_etag[e_tag]
            return self._joints[jt_id]
        else:
            return None