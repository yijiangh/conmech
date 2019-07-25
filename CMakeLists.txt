cmake_minimum_required(VERSION 3.1)
project(conmech)

# All dependencies that are downloaded as cmake projects and tested on the auto-builds are ON
# (by default, all build options are off)
option(CONMECH_WITH_PYTHON    "Use Python"    ON)
### End

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

include(conmech)
