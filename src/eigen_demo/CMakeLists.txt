add_definitions("-std=c++11")

include_directories(${EIGEN_INCLUDE_DIRS})

add_executable(ed_test test/eigen_demo_test.cpp)

#add_library(eigen_demo
#
#	)
#
#target_link_libraries(sc_test stiffness_checker)
