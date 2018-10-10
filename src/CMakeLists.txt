include_directories(
  ${CMAKE_SOURCE_DIR}/ext
  ${CMAKE_SOURCE_DIR}/src
)

add_subdirectory(eigen_demo)
#add_subdirectory(stiffness_checker)

if(BUILD_PYBINDINGS)
	message(STATUS "Building python bindings.")
	add_subdirectory(bindings/pyconmech)
endif(BUILD_PYBINDINGS)