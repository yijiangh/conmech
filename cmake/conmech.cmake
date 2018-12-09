cmake_minimum_required(VERSION 3.1)

# https://github.com/libigl/libigl/issues/751
# http://lists.llvm.org/pipermail/llvm-commits/Week-of-Mon-20160425/351643.html
if(APPLE)
    if(NOT CMAKE_LIBTOOL)
        find_program(CMAKE_LIBTOOL NAMES libtool)
    endif()
    if(CMAKE_LIBTOOL)
        set(CMAKE_LIBTOOL ${CMAKE_LIBTOOL} CACHE PATH "libtool executable")
        message(STATUS "Found libtool - ${CMAKE_LIBTOOL}")
        get_property(languages GLOBAL PROPERTY ENABLED_LANGUAGES)
        foreach(lang ${languages})
            # Added -c
            set(CMAKE_${lang}_CREATE_STATIC_LIBRARY
                    "${CMAKE_LIBTOOL} -c -static -o <TARGET> <LINK_FLAGS> <OBJECTS> ")
        endforeach()
    endif()
endif()

### Available options ###
option(CONMECH_WITH_OPENGL            "Use OpenGL"                   OFF)
option(CONMECH_WITH_OPENGL_GLFW       "Use GLFW"                     OFF)
option(CONMECH_WITH_OPENGL_GLFW_IMGUI "Use ImGui"                    OFF)
option(CONMECH_WITH_PYTHON            "Use Python"                   OFF)

################################################################################

### Configuration
set(CONMECH_ROOT "${CMAKE_CURRENT_LIST_DIR}/..")
set(CONMECH_SOURCE_DIR "${CONMECH_ROOT}/src")
set(CONMECH_EXTERNAL "${CONMECH_ROOT}/ext")

# Download and update 3rdparty libraries
#list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})
#include(ConmechDownloadExternal)

find_package(Eigen REQUIRED)

set(RAPIDJSON_INCLUDE_DIRS ${CONMECH_EXTERNAL}/rapidjson/include)

add_subdirectory(${CONMECH_EXTERNAL})
add_subdirectory(${CONMECH_SOURCE_DIR})

################################################################################