#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/iostream.h>
#include "pybind11_json/pybind11_json.hpp"
namespace py = pybind11;

#include <stiffness_checker/Stiffness.h>
#include <stiffness_checker/StiffnessIO.h>

PYBIND11_MODULE(_pystiffness_checker, m)
{
  m.doc() = "conmech stiffness checker";

  // TODO: https://pybind11.readthedocs.io/en/stable/advanced/classes.html#custom-constructors
  // https://github.com/jpanetta/MeshFEM/blob/master/src/python_bindings/MeshFactory.hh

  using Stiffness = conmech::stiffness_checker::Stiffness;
  py::class_<Stiffness>(m, "_StiffnessChecker", "stiffness checker pybind")
    // Bind the factory function as a constructor:
    // .def(py::init(&Stiffness::create))

    .def(py::init<const std::string&, const bool&, const std::string&, const bool&>(),
         py::arg("json_file_path"), py::arg("verbose") = false,
         py::arg("model_type") = "frame", py::arg("output_json") = false)

    /**
     * @brief Bind a lambda function returning a pointer wrapped in a holder:
     *  needed for converting py::object to nl::json
     * 
     */
    .def(py::init([](const Eigen::MatrixXd& V, const Eigen::MatrixXi& E, const Eigen::MatrixXi& Fixities, 
                     const std::vector<py::object>& mat_dicts, const bool& verbose, const std::string& model_type, const bool& output_json)
      {
        namespace nl = nlohmann;
        using Material = conmech::material::Material;
        // https://github.com/pybind/pybind11_json
        // nl::json j = obj;  // Automatic py::object->nl::json conversion
        // py::object o = j;  // Automatic nl::json->py::object conversion
        std::vector<Material> mats;
        for (py::object mat_pyobj : mat_dicts)
        {
          nl::json json_mat = mat_pyobj;
          mats.push_back(Material(json_mat));
        }
        return std::unique_ptr<Stiffness>(new Stiffness(V, E, Fixities, mats, verbose, model_type, output_json));
      }),// end lambda function
      py::arg("vertices"), py::arg("elements"), py::arg("fixities"), py::arg("material_dicts"),
      py::arg("verbose") = false, py::arg("model_type") = "frame", py::arg("output_json") = false)

    .def("set_self_weight_load", [](conmech::stiffness_checker::Stiffness &cm, 
                                    const bool &include_self_weight = true, 
                                    const std::vector<double> &gravity_direction = std::vector<double>{0,0,-1.0})
    {
      cm.setSelfWeightNodalLoad(include_self_weight);
      Eigen::VectorXd gravity_direction_vec = Eigen::VectorXd(int(gravity_direction.size()));
      for (int i=0; i<gravity_direction.size(); i++) {
        gravity_direction_vec[i] = gravity_direction[i];
      }
      cm.setGravityDirection(gravity_direction_vec);
    },
    py::arg("include_self_weight") = true, py::arg("gravity_direction") = std::vector<double>{0,0,-1.0})

    // input: #nL x 7 numpy matrix:
    // e.g.
    // L = np.zeros(nL, 7)
    // L[i,0] = node_id
    // L[i,1:7] = [0,0,-1,0,0,0]
    .def("set_load", &conmech::stiffness_checker::Stiffness::setLoad, py::arg("nodal_forces"), "input: #nL x 7 numpy matrix")

    .def("set_uniformly_distributed_loads", &conmech::stiffness_checker::Stiffness::setUniformlyDistributedLoad, py::arg("element_load_density"))

    // note: file_path must be attached with a path separator
    .def("set_output_json_path", &conmech::stiffness_checker::Stiffness::setOutputJsonPath,
      py::arg("file_path"), py::arg("file_name"))

    // output_json = True / False
    .def("set_output_json", &conmech::stiffness_checker::Stiffness::setOutputJson,
      py::arg("output_json") = false)

    // trans unit: meter, rot unit: rad
    .def("set_nodal_displacement_tol", &conmech::stiffness_checker::Stiffness::setNodalDisplacementTolerance,
      py::arg("transl_tol"), py::arg("rot_tol"))

    .def("solve",
      [](conmech::stiffness_checker::Stiffness &cm, const std::vector<int> &exist_element_ids = std::vector<int>(), const bool &if_cond_num = true)
      {
        if (exist_element_ids.empty()) 
        {
          return cm.solve(if_cond_num);
        } else 
        {
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
          for (int i=0;i<cm.nE();i++) tmp_existing_ids.push_back(i);
        }

        Eigen::VectorXd tot_pt_load;
        cm.getExternalNodalLoad(tot_pt_load);

        Eigen::VectorXd element_lumped_load;
        cm.getUniformlyDistributedLumpedLoad(tmp_existing_ids, element_lumped_load);
        tot_pt_load += element_lumped_load;

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
        for (int i=0;i<cm.nE();i++) tmp_existing_ids.push_back(i);
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
        return std::make_tuple(cm.nV(), cm.nE());
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

