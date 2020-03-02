cmake_minimum_required(VERSION 3.1)

################################################################################

### Configuration
set(CONMECH_ROOT "${CMAKE_CURRENT_LIST_DIR}/..")
set(CONMECH_SOURCE_DIR "${CONMECH_ROOT}/src")
set(CONMECH_EXTERNAL "${CONMECH_ROOT}/ext")

# Download and update 3rdparty libraries
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})

# C++11/14 features
include(CXXFeatures)

# https://github.com/jpanetta/MeshFEM/blob/master/CMakeLists.txt#L67
# We need -fPIC when compiling our libraries and our dependencies for
# the python bindings (shared libraries) to link.
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

################################################################################
# Dependencies
################################################################################

# Download and define targets for third-party dependencies
include(ConmechDependencies)

################################################################################
# Subdirectories
################################################################################

add_subdirectory(${CONMECH_SOURCE_DIR})

################################################################################
# Unit tests
################################################################################

if(CONMECH_BUILD_TESTS)
	include(CTest)
	enable_testing()
	add_subdirectory(tests/cpp)
endif()