cmake_minimum_required(VERSION 3.1)
project(conmech)

# TODO: All dependencies that are downloaded as cmake projects and tested on the auto-builds are ON
# (by default, all build options are off)
# https://github.com/libigl/libigl/blob/f6b406427400ed7ddb56cfc2577b6af571827c8c/cmake/LibiglDownloadExternal.cmake
# https://github.com/jpanetta/MeshFEM/tree/master/cmake
option(CONMECH_WITH_PYTHON    "Use Python"    OFF)
option(CONMECH_BUILD_TESTS    "Build tests"   ON)
### End

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

include(conmech)

################################################################################
# Unit tests
################################################################################

if(CONMECH_BUILD_TESTS)
	include(CTest)
	enable_testing()
	add_subdirectory(tests/cpp)
endif()

# if(LIBIGL_WITH_PYTHON)
# 	add_subdirectory(python)
# endif()