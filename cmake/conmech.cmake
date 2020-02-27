cmake_minimum_required(VERSION 3.1)

################################################################################

### Configuration
set(CONMECH_ROOT "${CMAKE_CURRENT_LIST_DIR}/..")
set(CONMECH_SOURCE_DIR "${CONMECH_ROOT}/src")
set(CONMECH_EXTERNAL "${CONMECH_ROOT}/ext")

# Download and update 3rdparty libraries
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})

################################################################################
# Dependencies
################################################################################

# Download and define targets for third-party dependencies
include(ConmechDependencies)

# for now, eigen is shipped with conmech
# find_package(Eigen REQUIRED)
# set(EIGEN_INCLUDE_DIRS ${CONMECH_EXTERNAL}/eigen)
# set(RAPIDJSON_INCLUDE_DIRS ${CONMECH_EXTERNAL}/rapidjson/include)

# add_subdirectory(${CONMECH_EXTERNAL})

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