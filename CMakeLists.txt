cmake_minimum_required(VERSION 2.8.12)
project(conmech)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules")
find_package(Eigen REQUIRED)

set(RAPIDJSON_INCLUDE_DIRS ${CMAKE_SOURCE_DIR}/ext/rapidjson/include)

#set(BUILD_PYBINDINGS FALSE)
set(BUILD_PYBINDINGS TRUE)

add_subdirectory(ext)
add_subdirectory(src)
