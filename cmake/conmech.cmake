cmake_minimum_required(VERSION 3.1)

### Configuration
set(CONMECH_ROOT "${CMAKE_CURRENT_LIST_DIR}/..")
set(CONMECH_SOURCE_DIR "${CONMECH_ROOT}/src")
set(CONMECH_EXTERNAL "${CONMECH_ROOT}/ext")

# Download and update 3rdparty libraries
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})
# include(ConmechDownloadExternal)

find_package(Eigen REQUIRED)
set(RAPIDJSON_INCLUDE_DIRS ${CONMECH_EXTERNAL}/rapidjson/include)

add_subdirectory(${CONMECH_EXTERNAL})
add_subdirectory(${CONMECH_SOURCE_DIR})

################################################################################
