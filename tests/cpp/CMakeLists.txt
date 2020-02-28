################################################################################
# Independent binaries
################################################################################
# https://github.com/jpanetta/MeshFEM/blob/master/tests/CMakeLists.txt

add_executable(unit_tests
	main.cpp
    test_memory_leak.cpp
    test_io.cpp
)

target_link_libraries(unit_tests PUBLIC
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
catch_discover_tests(unit_tests)

add_test(unit_tests ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/unit_tests)