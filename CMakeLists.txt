cmake_minimum_required(VERSION 3.1)
project(conmech)

# # Check if MeshFEM is being used directly or via add_subdirectory
# set(MeshFEM_MASTER_PROJECT OFF)
# if (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
#   set(MeshFEM_MASTER_PROJECT ON)
# endif()

option(CONMECH_WITH_PYTHON    "Use Python"    OFF)
option(CONMECH_BUILD_TESTS    "Build tests"   ON)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

include(conmech)