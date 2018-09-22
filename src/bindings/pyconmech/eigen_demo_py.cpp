#include <pybind11/pybind11.h>

#include "eigen_demo/eigen_demo.cpp"

namespace py = pybind11;

namespace conmech
{

namespace pyconmech
{

PYBIND11_MODULE(eigen_demo_py, m
)
{
py::class_<Pet>(m,
"Pet")
.
def (py::init<const std::string &>())
.def("setName", &Pet::setName)
.def("getName", &Pet::getName)

.def_property("name", &Pet::getName, &Pet::setName)
.def_readwrite("file_path", &Pet::file_path)

.def("__repr__",
[](
const Pet &a
)
{
return "<example.Pet named '" + a.
getName()
+ "'>";
}
);
}

} // namespace pyconmech
} // namespace conmech
