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
.def(py::init<const std::string&, bool>(), py::arg("json_file_path"), py::arg("verbose") = false)

.def("set_self_weight_nodal_load", &conmech::stiffness_checker::Stiffness::setSelfWeightNodalLoad)

.def("set_nodal_load", &conmech::stiffness_checker::Stiffness::setNodalLoad,
py::arg("nodal_forces"), py::arg("include_self_weight") = false)

.def("solve", py::overload_cast<const std::vector<int> &, const bool&>(&conmech::stiffness_checker::Stiffness::solve),
py::arg("exist_element_ids"), py::arg("if_cond_num") = true)

.def("solve", py::overload_cast<const bool&>(&conmech::stiffness_checker::Stiffness::solve),
    py::arg("if_cond_num") = true)
    ;

}

} // namespace pyconmech
} // namespace conmech