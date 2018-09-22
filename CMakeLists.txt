cmake_minimum_required(VERSION 2.8.12)
project(conmech)

#list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules")

#find_package(Eigen REQUIRED)

#set(JSON_INCLUDE_DIRS ${CMAKE_SOURCE_DIR}/ext/json/include)

add_subdirectory(ext)
add_subdirectory(src)
