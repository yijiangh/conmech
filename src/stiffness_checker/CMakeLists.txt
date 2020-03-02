################################################################################

add_library(conmech_stiffness_checker
	Stiffness.cpp
	StiffnessSolver.cpp
	StiffnessIO.cpp
	StiffnessShape.cpp
	Util.cpp
	Timer.cpp
	Material.cpp
	)
# # Export conmech_stiffness_checker as conmech::stiffness_checker
# set_property(TARGET conmech_stiffness_checker PROPERTY EXPORT_NAME conmech::stiffness_checker)

# Include directories
target_include_directories(conmech_stiffness_checker PUBLIC ..)

# Use C++11
target_compile_features(conmech_stiffness_checker PUBLIC ${CXX11_FEATURES})

# header-only libraries linking
# https://github.com/jpanetta/MeshFEM/blob/master/src/lib/MeshFEM/CMakeLists.txt
target_link_libraries(conmech_stiffness_checker
    PUBLIC
        Eigen3::Eigen
		json::json
    # PRIVATE
    #     warnings::all
)
add_library(conmech::stiffness_checker ALIAS conmech_stiffness_checker)