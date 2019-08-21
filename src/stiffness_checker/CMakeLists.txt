# add_definitions("-std=c++11")
# TODO: lapack is needed for condition number checking
#find_package(LAPACK REQUIRED)

include_directories(${EIGEN_INCLUDE_DIRS} ${RAPIDJSON_INCLUDE_DIRS})

# add_executable(sc_test test/stiffness_checker_test.cpp)

add_library(stiffness_checker_lib
	Frame.cpp
	Util.cpp
	Stiffness.cpp
	StiffnessSolver.cpp
	Timer.cpp
	StiffnessIO.cpp
	)

# auto deducing CXX features
# https://stackoverflow.com/questions/10851247/how-to-activate-c-11-in-cmake
# https://stackoverflow.com/questions/10984442/how-to-detect-c11-support-of-a-compiler-with-cmake/20165220#20165220
target_compile_features(stiffness_checker_lib PRIVATE cxx_range_for)

# target_link_libraries(sc_test stiffness_checker_lib)

if(CONMECH_WITH_PYTHON)
	message(STATUS "Building stiffness_checker's python bindings.")
	add_subdirectory(bindings)
endif(CONMECH_WITH_PYTHON)