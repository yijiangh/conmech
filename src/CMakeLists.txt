add_subdirectory(stiffness_checker)

################################################################################
# Python Bindings
################################################################################

if(CONMECH_WITH_PYTHON)
	add_subdirectory(pyconmech_bindings)
endif()