################################################################################
# Independent binaries
################################################################################
# https://github.com/jpanetta/MeshFEM/blob/master/tests/CMakeLists.txt
# https://github.com/libigl/libigl/blob/master/tests/CMakeLists.txt

# Create test executable
add_executable(conmech_tests
	main.cpp
	test_common.h
    test_memory_leak.cpp
    test_io.cpp
)

# Set DATA_DIR definition
set(CONMECH_TEST_DATA "${CMAKE_CURRENT_SOURCE_DIR}/../test_data/")
target_compile_definitions(conmech_tests PUBLIC -DCONMECH_DATA_DIR="${CONMECH_TEST_DATA}")

target_link_libraries(conmech_tests PUBLIC
	Catch2::Catch2
	json::json
	stiffness_checker_lib
)

# if (MSVC)
# else ()
# target_compile_options(unit_tests PRIVATE -Wno-unused-parameter)
# endif()

# Register tests
include(CTest)
include(Catch)
# set(PARSE_CATCH_TESTS_ADD_TO_CONFIGURE_DEPENDS ON)
catch_discover_tests(conmech_tests)

add_test(conmech_tests ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/conmech_tests)