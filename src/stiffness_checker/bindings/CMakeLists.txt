# needed for pybind11 to find eigen path
include_directories(${EIGEN_INCLUDE_DIRS})

pybind11_add_module(_pystiffness_checker
	pybind_stiffness_checker.cpp
	)

target_link_libraries(_pystiffness_checker PRIVATE stiffness_checker_lib)

# this is needed for linking a static lib into a dynamic lib
# see:
# https://stackoverflow.com/questions/38296756/what-is-the-idiomatic-way-in-cmake-to-add-the-fpic-compiler-option
# https://github.com/pybind/pybind11/issues/466
set_target_properties(stiffness_checker_lib PROPERTIES POSITION_INDEPENDENT_CODE TRUE)