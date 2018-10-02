#include <pybind11/pybind11.h>

#include "stiffness_checker/Stiffness.hpp"

namespace py = pybind11;

namespace conmech
{
namespace pyconmech
{

PYBIND11_MODULE(conmech_py, m)
{

py::class_<StiffnessChecker>(m,"stiffness_checker")
.def(py::init<std::string, bool>(), "verbose"_a=false)
.def("check_deformation",
py::overload_cast<const std::vector<int> &>(&conmech::stiffness_checker::StiffnessChecker::checkDeformation))
.def("check_deformation",
py::overload_cast<const std::vector<int> &, std::vector<double> &>(&conmech::stiffness_checker::StiffnessChecker::checkDeformation));

}

} // namespace pyconmech
} // namespace conmech