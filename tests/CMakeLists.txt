#
# Build examples.
#
# Requires CONMECH_BUILD_TESTS to be defined 'true', see ../CMakeLists.txt.
#

# cmake_minimum_required( VERSION 3.0 )
# project( conmech_tests CXX )

message( STATUS "Tests included" )

add_executable(conmech_tests
hello_catch.cpp
)
# target_link_libraries(inja_test PRIVATE inja)

add_test(conmech_tests ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/conmech_tests)

# add_executable(inja_test
# test/unit.cpp
# test/unit-files.cpp
# test/unit-renderer.cpp
# )
# target_link_libraries(inja_test PRIVATE inja)

# add_test(inja_test ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/inja_test)


# add_library(single_inja INTERFACE)
# target_compile_features(single_inja INTERFACE cxx_std_11)
# target_include_directories(single_inja INTERFACE single_include include third_party/include)

# add_executable(single_inja_test
# test/unit.cpp
# test/unit-files.cpp
# test/unit-renderer.cpp
# )
# target_link_libraries(single_inja_test PRIVATE single_inja)

# add_test(single_inja_test ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/single_inja_test)