
pybind11_add_module(conmech_py
	stiffness_checker/stiffness_checker_py.cpp)

#target_link_libraries(eigen_demo_py PRIVATE eigen_demo)

# this is needed for linking a static lib into a dynamic lib
# see:
# https://stackoverflow.com/questions/38296756/what-is-the-idiomatic-way-in-cmake-to-add-the-fpic-compiler-option
# https://github.com/pybind/pybind11/issues/466
#set_target_properties(stiffness_checker PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
