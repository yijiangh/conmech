import numpy as np
from .stiffness_base import StiffnessBase

class NumpyStiffness(StiffnessBase):
    def __init__(self, vertices, elements, fixities, material_dicts, 
        verbose=False, model_type='frame', output_json=False):
        super(NumpyStiffness, self).__init__(vertices, elements, fixities, material_dicts, \
            verbose, model_type, output_json)
        # default displacement tolerance
        self._trans_tol = 1e-3
        self._rot_tol = np.pi/180 * 3

    @classmethod
    def from_json(cls, json_file_path=None, verbose=False):
        raise NotImplementedError()

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