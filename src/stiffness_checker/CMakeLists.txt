################################################################################

add_library(stiffness_checker_lib
	Stiffness.cpp
	StiffnessSolver.cpp
	StiffnessIO.cpp
	Util.cpp
	Timer.cpp
	Material.cpp
	)

# Include directories
target_include_directories(stiffness_checker_lib PUBLIC ..)

# Use C++11
# target_compile_features(${PROJECT_NAME} PUBLIC ${CXX11_FEATURES})

# header-only libraries linking
# https://github.com/jpanetta/MeshFEM/blob/master/src/lib/MeshFEM/CMakeLists.txt
target_link_libraries(stiffness_checker_lib
    PUBLIC
        Eigen3::Eigen
		json::json
    # PRIVATE
    #     warnings::all
)

if(CONMECH_WITH_PYTHON)
	message(STATUS "Building stiffness_checker's python bindings.")
	add_subdirectory(bindings)
endif()