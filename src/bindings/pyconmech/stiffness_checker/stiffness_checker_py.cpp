#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include <stiffness_checker/Stiffness.h>

namespace py = pybind11;

namespace conmech
{
namespace pyconmech
{

PYBIND11_MODULE(conmech_py, m)
{

py::class_<conmech::stiffness_checker::Stiffness>(m,"stiffness_checker")
.def(py::init<const std::string&, bool, const std::string&>(),
  py::arg("json_file_path"), py::arg("verbose") = false, py::arg("model_type") = "frame")

.def("set_self_weight_load", &conmech::stiffness_checker::Stiffness::setSelfWeightNodalLoad,
  py::arg("include_sw"))

.def("set_load", &conmech::stiffness_checker::Stiffness::setLoad, py::arg("nodal_forces"))

.def("set_nodal_displacement_tol", &conmech::stiffness_checker::Stiffness::setNodalDisplacementTolerance,
  py::arg("transl_tol"), py::arg("rot_tol"))

.def("solve", py::overload_cast<const std::vector<int> &, const bool&>(&conmech::stiffness_checker::Stiffness::solve),
py::arg("exist_element_ids"), py::arg("if_cond_num") = true)

// directly computing the whole structure
.def("solve", py::overload_cast<const bool&>(&conmech::stiffness_checker::Stiffness::solve),
    py::arg("if_cond_num") = true)
    ;

}

} // namespace pyconmech
} // namespace conmech