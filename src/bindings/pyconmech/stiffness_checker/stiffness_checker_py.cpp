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
.def(py::init<std::string>())
.def("check_deformation", &conmech::stiffness_checker::StiffnessChecker::checkDeformation);

}

} // namespace pyconmech
} // namespace conmech