#include <pybind11/pybind11.h>

#include "eigen_demo/eigen_demo.cpp"

namespace py = pybind11;

namespace conmech
{

namespace pyconmech
{

PYBIND11_MODULE(eigen_demo_py, m)
{

py::class_<EigenSolveDemo>(m,"EigenSolveDemo")
.def(py::init<>())
.def("testEigen", &EigenSolveDemo::testEigen);

}

} // namespace pyconmech
} // namespace conmech
