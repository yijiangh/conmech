add_definitions("-std=c++11")

include_directories(${EIGEN_INCLUDE_DIRS})

add_executable(eigen_demo_test eigen_demo_test.cpp)

add_library(eigen_demo
		eigen_demo.cpp
	)

target_link_libraries(eigen_demo_test eigen_demo)
