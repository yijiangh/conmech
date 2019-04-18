include_directories(
		${CONMECH_EXTERNAL}
		${CONMECH_SOURCE_DIR})

add_subdirectory(eigen_demo)
add_subdirectory(stiffness_checker)

if(CONMECH_WITH_PYTHON)
	message(STATUS "Building python bindings.")
	add_subdirectory(bindings/pyconmech)
endif(CONMECH_WITH_PYTHON)
