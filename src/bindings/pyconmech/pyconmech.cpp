#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/iostream.h>

#include <stiffness_checker/Stiffness.h>
#include <stiffness_checker/StiffnessIO.h>
#include <eigen_demo/eigen_demo.h>

namespace py = pybind11;

namespace conmech
{
namespace pyconmech
{

PYBIND11_MODULE(pyconmech, m)
{
    py::class_<conmech::stiffness_checker::Stiffness>(m,"stiffness_checker")
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

    .def("solve", py::overload_cast<const std::vector<int> &, const bool&>(&conmech::stiffness_checker::Stiffness::solve),
    py::arg("exist_element_ids"), py::arg("if_cond_num") = true)

    // directly computing the whole structure
    .def("solve", py::overload_cast<const bool&>(&conmech::stiffness_checker::Stiffness::solve),
    py::arg("if_cond_num") = true)

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
        if(cm.hasStoredResults()) {
            cm.getMaxNodalDeformation(max_trans, max_rot);
        }
        return std::make_tuple(max_trans, max_rot);
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

m.def("parse_load_case_from_json", [](std::string file_path)
{
    Eigen::MatrixXd Load;
    bool include_sw;
    conmech::stiffness_checker::parseLoadCaseJson(file_path, Load, include_sw);
    return std::make_tuple(Load, include_sw);
});

py::class_<EigenSolveDemo>(m,"eigen_solve_demo")
  .def(py::init<>())
  .def("test_eigen", &EigenSolveDemo::testEigen)
  ; // end eigen solve demo

// TODO: redirecting C++ streams
// https://pybind11.readthedocs.io/en/stable/reference.html#redirecting-c-streams
py::add_ostream_redirect(m, "ostream_redirect");

} // end pyconmech def

} // namespace pyconmech
} // namespace conmech
