add_definitions("-std=c++11")

#find_package(LAPACK REQUIRED)

include_directories(${EIGEN_INCLUDE_DIRS} ${RAPIDJSON_INCLUDE_DIRS})

#add_executable(sc_test test/stiffness_checker_test.cpp)

add_library(stiffness_checker
#	StiffnessChecker.cpp
#	Stiffness.cpp
#	StiffnessSolver.cpp
#	StiffnessIO.cpp
#	DualGraph.cpp
	Frame.cpp
#	CoordTrans.cpp
#	Timer.cpp
	)

#target_link_libraries(sc_test stiffness_checker)