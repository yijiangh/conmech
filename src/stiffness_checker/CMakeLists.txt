################################################################################

add_library(stiffness_checker_lib
	Stiffness.cpp
	StiffnessSolver.cpp
	StiffnessIO.cpp
	Util.cpp
	Timer.cpp
	)

# Include directories
# include_directories(${EIGEN_INCLUDE_DIRS} ${RAPIDJSON_INCLUDE_DIRS})
target_include_directories(stiffness_checker_lib 
	PUBLIC 
	..
	)

# # auto deducing CXX features
# # https://stackoverflow.com/questions/10851247/how-to-activate-c-11-in-cmake
# # https://stackoverflow.com/questions/10984442/how-to-detect-c11-support-of-a-compiler-with-cmake/20165220#20165220
# target_compile_features(stiffness_checker_lib PRIVATE cxx_range_for)

# # Use C++11
# target_compile_features(${PROJECT_NAME} PUBLIC ${CXX11_FEATURES})

# header-only libraries linking
# https://github.com/jpanetta/MeshFEM/blob/master/src/lib/MeshFEM/CMakeLists.txt
target_link_libraries(stiffness_checker_lib
    PUBLIC
        Eigen3::Eigen
		json::json
		rapidjson::rapidjson
    # PRIVATE
    #     warnings::all
)

if(CONMECH_WITH_PYTHON)
	message(STATUS "Building stiffness_checker's python bindings.")
	add_subdirectory(bindings)
endif()