#add_subdirectory(rapidjson)

if(BUILD_PYBINDINGS)
	add_subdirectory(pybind11)
endif(BUILD_PYBINDINGS)