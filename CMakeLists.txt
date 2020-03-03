cmake_minimum_required(VERSION 3.1)
project(conmech)

option(CONMECH_WITH_PYTHON    "Use Python"    ON)
option(CONMECH_BUILD_TESTS    "Build tests"   ON)

### conditionally compile certain modules depending on libraries found on the system
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

include(conmech)