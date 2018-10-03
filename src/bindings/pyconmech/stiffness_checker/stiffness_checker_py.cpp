#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "stiffness_checker/StiffnessChecker.hpp"

namespace py = pybind11;

namespace conmech
{
namespace pyconmech
{

PYBIND11_MODULE(conmech_py, m)
{

py::class_<conmech::stiffness_checker::StiffnessChecker>(m,"stiffness_checker")
.def(py::init<std::string, bool>())
.def("check_deformation",
py::overload_cast<const std::vector<int> &>(&conmech::stiffness_checker::StiffnessChecker::checkDeformation))
.def("check_deformation",
py::overload_cast<const std::vector<int> &, std::vector<double> &>
    (&conmech::stiffness_checker::StiffnessChecker::checkDeformation), py::return_value_policy::copy);

}

} // namespace pyconmech
} // namespace conmech