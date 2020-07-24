import os
from pyconmech.frame_analysis.frame_file_io import read_frame_json, read_load_case_json
    
##################################################
# Stiffness Backend Template

class StiffnessBase(object):
    def __init__(self, nodes, elements, fixities, material_dicts, 
        verbose=False, model_type='frame', output_json=False):
        self._nodes = vertices
        self._elements = elements
        self._fixties = fixities
        self._material_dicts = material_dicts
        self._verbose = verbose
        self._model_type = model_type
        self._output_json = output_json
        self._trans_tol = None
        self._rot_tol = None

        if len(self._fixties) == 0:
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
        node_points, elements, fix_specs, model_type, material_dicts, _, unit = \
            read_frame_json(json_file_path, verbose=verbose)
        return cls(node_points, elements, fix_specs, material_dicts, 
            model_type=model_type, verbose=verbose)

    #######################
    # Load input

    def set_self_weight_load(self, include_self_weight=True, gravity_direction=[0,0,-1]):
        raise NotImplementedError()

    def set_load(self, nodal_forces):
        """input: #nL x 7 numpy matrix
        """
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