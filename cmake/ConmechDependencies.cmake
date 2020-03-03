# Prepare dependencies
#
# For each third-party library, if the appropriate target doesn't exist yet,
# download it via external project, and add_subdirectory to build it alongside
# this project.

# Download and update 3rdparty libraries
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR})
list(REMOVE_DUPLICATES CMAKE_MODULE_PATH)
include(ConmechDownloadExternal)

################################################################################
# Required libraries
################################################################################
# https://cmake.org/cmake/help/latest/manual/cmake-buildsystem.7.html#interface-libraries

# Eigen3 library
if(TARGET Eigen3::Eigen)
    # If an imported target already exists, use it
    target_link_libraries(conmech_eigen INTERFACE Eigen3::Eigen)
else()
    add_library(conmech_eigen INTERFACE)
    conmech_download_eigen()
#   target_include_directories(igl_common SYSTEM INTERFACE
#     $<BUILD_INTERFACE:${LIBIGL_EXTERNAL}/eigen>
#     $<INSTALL_INTERFACE:include>
#   )
    target_include_directories(conmech_eigen SYSTEM INTERFACE ${CONMECH_EXTERNAL}/eigen)
    add_library(Eigen3::Eigen ALIAS conmech_eigen)
endif()

# json library
if(NOT TARGET json::json)
    add_library(conmech_json INTERFACE)
    conmech_download_json()
    target_include_directories(conmech_json SYSTEM INTERFACE ${CONMECH_EXTERNAL}/json/include)
    add_library(json::json ALIAS conmech_json)
endif()

# Catch2
if(NOT TARGET Catch2::Catch2 AND (CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR))
    conmech_download_catch()
    add_subdirectory(${CONMECH_EXTERNAL}/Catch2)
    list(APPEND CMAKE_MODULE_PATH ${CONMECH_EXTERNAL}/Catch2/contrib)
endif()

################################################################################
### Download the python part ###
if(CONMECH_WITH_PYTHON)
    if(NOT TARGET pybind11)
        conmech_download_project(pybind11
            GIT_REPOSITORY https://github.com/pybind/pybind11.git
            GIT_TAG 97784dad3e518ccb415d5db57ff9b933495d9024
        )
    endif()
    add_subdirectory(${CONMECH_EXTERNAL}/pybind11)
endif()