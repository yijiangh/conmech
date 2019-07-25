if(CONMECH_WITH_PYTHON)
	# pybind options
	# https://pybind11.readthedocs.io/en/stable/compiling.html

	# check avaible C++ features
	# message("Your C++ compiler supports these C++ features:")
	# foreach(i ${CMAKE_CXX_COMPILE_FEATURES})
	#   message("${i}")
	# endforeach()

	add_subdirectory(pybind11)
endif(CONMECH_WITH_PYTHON)
