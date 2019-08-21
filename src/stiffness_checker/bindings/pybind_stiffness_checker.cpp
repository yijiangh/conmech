#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/iostream.h>

#include <stiffness_checker/Stiffness.h>
#include <stiffness_checker/StiffnessIO.h>

namespace py = pybind11;

namespace conmech {
namespace pyconmech {

PYBIND11_MODULE(_pystiffness_checker, m)
{
    py::class_<conmech::stiffness_checker::Stiffness>(m, "_stiffness_checker")
    .def(py::init<const std::string&, bool, const std::string&, bool>(),
    py::arg("json_file_path"), py::arg("verbose") = false,
    py::arg("model_type") = "frame", py::arg("output_json") = false)

    .def("set_self_weight_load", &conmech::stiffness_checker::Stiffness::setSelfWeightNodalLoad,
    py::arg("include_sw") = false)

    // input: 1 x 7 numpy matrix:
    // e.g.
    // L = np.zeros(nL, 7)
    // L[i,0] = node_id
    // L[i,1:7] = [0,0,-1,0,0,0]
    .def("set_load", &conmech::stiffness_checker::Stiffness::setLoad, py::arg("nodal_forces"))

    // note: file_path must be attached with a path separator
    .def("set_output_json_path", &conmech::stiffness_checker::Stiffness::setOutputJsonPath,
    py::arg("file_path"), py::arg("file_name"))

    // output_json = True / False
    .def("set_output_json", &conmech::stiffness_checker::Stiffness::setOutputJson,
    py::arg("output_json") = false)

    // trans unit: meter, rot unit: rad
    .def("set_nodal_displacement_tol", &conmech::stiffness_checker::Stiffness::setNodalDisplacementTolerance,
    py::arg("transl_tol"), py::arg("rot_tol"))


    // element length L: m
    // E, G: kN/m^2
    // Force: kN
    // Moment: kN m
    .def("solve",
    [](conmech::stiffness_checker::Stiffness &cm, const std::vector<int> &exist_element_ids = std::vector<int>(), const bool &if_cond_num = true)
    {
      // TODO: sanity check existing_ids within range
      if (exist_element_ids.empty()) {
        return cm.solve(if_cond_num);
      } else {
        return cm.solve(exist_element_ids, if_cond_num);
      }
    },
    py::arg("exist_element_ids") = std::vector<int>(),
    py::arg("if_cond_num") = true)

    // for my own education:
    // C++ lambda function stateful/stateless enclosure
    // https://arne-mertz.de/2015/10/new-c-features-lambdas/
    // https://arne-mertz.de/2015/11/lambdas-part-2-capture-lists-and-stateful-closures/
    // Get all lumped nodal load at the node 
    .def("get_lumped_nodal_loads",
    [](conmech::stiffness_checker::Stiffness &cm, const std::vector<int> &existing_ids = std::vector<int>())
    {
      std::vector<int> tmp_existing_ids = existing_ids;
      // TODO: sanity check existing_ids within range
      if (tmp_existing_ids.empty()) {
        for (int i=1;i<cm.getTotalNumOfElements();i++) tmp_existing_ids.push_back(i);
      }
      Eigen::VectorXd tot_pt_load;
      cm.getExternalNodalLoad(tot_pt_load);
      if (cm.isIncludeSelfWeightLoad()) {
        Eigen::VectorXd sw_nodal_loads;
        cm.getSelfWeightNodalLoad(tmp_existing_ids, sw_nodal_loads);
        tot_pt_load += sw_nodal_loads;
      }
      return tot_pt_load;
    },
    py::arg("existing_ids") = std::vector<int>())

    .def("get_gravity_nodal_loads",
    [](conmech::stiffness_checker::Stiffness &cm, const std::vector<int> &existing_ids = std::vector<int>())
    {
      std::vector<int> tmp_existing_ids = existing_ids;
      // TODO: sanity check existing_ids within range
      if (existing_ids.empty()) {
        for (int i=1;i<cm.getTotalNumOfElements();i++) tmp_existing_ids.push_back(i);
      }
      Eigen::VectorXd sw_nodal_loads;
      cm.getSelfWeightNodalLoad(tmp_existing_ids, sw_nodal_loads);
      return sw_nodal_loads;
    },
    py::arg("existing_ids") = std::vector<int>())

    // get a list of elemental stiffness matrix (R * K_{eL} * R^T)
    // all in global coordinate system, 12 x 12 matrix
    .def("get_element_stiffness_matrices",
    [](conmech::stiffness_checker::Stiffness &cm)
    {
      std::vector<Eigen::MatrixXd> element_stiffness_mats;
      cm.getElementStiffnessMatrices(element_stiffness_mats);
      return element_stiffness_mats;
    })

    // get a list of elemental local to global rotational matrix, 12 x 12 matrix
    .def("get_element_local2global_rot_matrices",
    [](conmech::stiffness_checker::Stiffness &cm)
    {
      std::vector<Eigen::MatrixXd> element_rot_mats;
      cm.getElementLocal2GlobalRotationMatrices(element_rot_mats);
      return element_rot_mats;
    })

    // returns a (N_element x (2*6)) map
    // element id -> dof id map
    .def("get_element2dof_id_map",
    [](conmech::stiffness_checker::Stiffness &cm)
    {
        return cm.getElement2DofIdMap();
    })

    // returns a (N_node x (6)) map
    // node id -> dof id map
    .def("get_node2dof_id_map",
    [](conmech::stiffness_checker::Stiffness &cm)
    {
        return cm.getNode2DofIdMap();
    })

    // check if have stored results
    .def("has_stored_result", &conmech::stiffness_checker::Stiffness::hasStoredResults)

    // get solved result
    .def("get_solved_results",
    // TODO: make this const
    [](conmech::stiffness_checker::Stiffness &cm)
    {
        Eigen::MatrixXd node_displ, fixities_reaction, element_reaction;
        bool pass_criteria = false;
        if(cm.hasStoredResults())
        {
            cm.getSolvedResults(node_displ, fixities_reaction, element_reaction, pass_criteria);
        }
        return std::make_tuple(pass_criteria, node_displ, fixities_reaction, element_reaction);
    })

    .def("get_max_nodal_deformation",
    [](conmech::stiffness_checker::Stiffness &cm)
    {
        double max_trans, max_rot;
        int max_trans_vid, max_rot_vid;
        if(cm.hasStoredResults()) {
            cm.getMaxNodalDeformation(max_trans, max_rot, max_trans_vid, max_rot_vid);
        }
        return std::make_tuple(max_trans, max_rot, max_trans_vid, max_rot_vid);
    })

    // return (n_node, n_elment) for the full structure
    .def("get_frame_stat",
    [](conmech::stiffness_checker::Stiffness &cm)
    {
        return std::make_tuple(cm.getTotalNumOfVertices(), cm.getTotalNumOfElements());
    })

    .def("get_compliance",
    [](conmech::stiffness_checker::Stiffness &cm)
    {
        double compliance;
        if(cm.hasStoredResults()) {
            cm.getSolvedCompliance(compliance);
        }
        return compliance;
    })

    .def("get_nodal_deformation_tol",
    [](conmech::stiffness_checker::Stiffness &cm)
    {
        return std::make_tuple(cm.getTransTol(), cm.getRotTol());
    })

    // return the original, undeformed beam
    .def("get_original_shape", &conmech::stiffness_checker::Stiffness::getOriginalShape,
    py::arg("disc") = 1, py::arg("draw_full_shape") = true)

    // return the cubic interpolated deform beam
    .def("get_deformed_shape", &conmech::stiffness_checker::Stiffness::getDeformedShape,
    py::arg("exagg_ratio") = 1.0, py::arg("disc") = 10)
    // py::return_value_policy::reference_internal)

  ; // end stiffness checker

  m.def("_parse_load_case_from_json", [](std::string file_path)
  {
      Eigen::MatrixXd Load;
      bool include_sw;
      conmech::stiffness_checker::parseLoadCaseJson(file_path, Load, include_sw);
      return std::make_tuple(Load, include_sw);
  });

  // TODO: redirecting C++ streams
  // https://pybind11.readthedocs.io/en/stable/reference.html#redirecting-c-streams
  py::add_ostream_redirect(m, "ostream_redirect");

} // end pyconmech def

} // namespace pyconmech
} // namespace conmech
