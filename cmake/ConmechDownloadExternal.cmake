################################################################################
# https://github.com/jpanetta/MeshFEM/blob/master/cmake/MeshFEMDownloadExternal.cmake
include(DownloadProject)

# With CMake 3.8 and above, we can hide warnings about git being in a
# detached head by passing an extra GIT_CONFIG option
if(NOT (${CMAKE_VERSION} VERSION_LESS "3.8.0"))
    set(CONMECH_EXTRA_OPTIONS "GIT_CONFIG advice.detachedHead=false")
else()
    set(CONMECH_EXTRA_OPTIONS "")
endif()

# Shortcut function
function(conmech_download_project name)
    download_project(
        PROJ         ${name}
        SOURCE_DIR   ${CONMECH_EXTERNAL}/${name}
        DOWNLOAD_DIR ${CONMECH_EXTERNAL}/.cache/${name}
        QUIET
        ${CONMECH_EXTRA_OPTIONS}
        ${ARGN}
    )
endfunction()

################################################################################

## Catch2
function(conmech_download_catch)
    conmech_download_project(Catch2
        URL     https://github.com/catchorg/Catch2/archive/v2.3.0.tar.gz
        URL_MD5 1fc90ff3b7b407b83057537f4136489e
    )
endfunction()

## Eigen
set(CONMECH_EIGEN_VERSION 3.3.7 CACHE STRING "Default version of Eigen used by conmech.")
function(conmech_download_eigen)
    conmech_download_project(eigen
		GIT_REPOSITORY https://github.com/eigenteam/eigen-git-mirror.git
		GIT_TAG        ${CONMECH_EIGEN_VERSION}
    )
endfunction()

## Json
function(conmech_download_json)
    conmech_download_project(json
        URL      https://github.com/nlohmann/json/releases/download/v3.7.3/include.zip
        URL_HASH SHA256=87b5884741427220d3a33df1363ae0e8b898099fbc59f1c451113f6732891014
    )
endfunction()

################################################################################

## Test data
# function(igl_download_test_data)
# 	igl_download_project_aux(test_data
# 		"${LIBIGL_EXTERNAL}/../tests/data"
# 		GIT_REPOSITORY https://github.com/libigl/libigl-tests-data
# 		GIT_TAG        adc66cabf712a0bd68ac182b4e7f8b5ba009c3dd
# 	)
# endfunction()