# -Ax64
# https://stackoverflow.com/questions/45397986/pybind11-python-is-64-bit-chosen-compiler-is-32-bito

pybind11_add_module(_pystiffness_checker
	pybind_stiffness_checker.cpp
	)

## Add dependencies
target_link_libraries(_pystiffness_checker PUBLIC conmech::stiffness_checker)

# set_target_properties(_pystiffness_checker PROPERTIES PREFIX "")
# set_target_properties(_pystiffness_checker PROPERTIES LIBRARY_OUTPUT_DIRECTORY $<1:${CONMECH_SOURCE_DIR}/pyconmech>)

# https://github.com/jpanetta/MeshFEM/blob/master/python/init_template.py
# set(PYTHON_DEPLOY_DIR "${CMAKE_SOURCE_DIR}/src/pyconmech/cpp_stiffness_checker")
# file(MAKE_DIRECTORY "${PYTHON_DEPLOY_DIR}")
# configure_file(${PROJECT_SOURCE_DIR}/src/pyconmech/init_template.py "${PYTHON_DEPLOY_DIR}/__init__.py")