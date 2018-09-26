include_directories(${EIGEN_INCLUDE_DIRS} ${RAPIDJSON_INCLUDE_DIRS})

add_executable(sc_test test/stiffness_checker_test.cpp)

add_library(stiffness_checker
	StiffnessIO.cpp
	DualGraph.cpp
	WireFrame.cpp
	)

target_link_libraries(sc_test stiffness_checker)