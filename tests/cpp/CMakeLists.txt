################################################################################
# Independent binaries
################################################################################

add_executable(unit_tests
	main.cpp
    hello_catch.cpp
    # memory_check.cpp
)

target_link_libraries(unit_tests PUBLIC
	# Catch2::Catch2
	stiffness_checker_lib
)

# if (MSVC)
# else ()
# target_compile_options(unit_tests PRIVATE -Wno-unused-parameter)
# endif()

# Register tests
include(CTest)
# include(Catch)
# catch_discover_tests(unit_tests)

add_test(unit_tests ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/unit_tests)

#################################################

# message( STATUS "Tests included" )

# add_executable(conmech_tests
# hello_catch.cpp
# )
# # target_link_libraries(inja_test PRIVATE inja)

# add_test(conmech_tests ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/conmech_tests)

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