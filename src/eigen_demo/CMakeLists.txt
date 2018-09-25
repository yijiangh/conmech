include_directories(${EIGEN_INCLUDE_DIRS} ${RAPIDJSON_INCLUDE_DIRS})

# build a C++ library
add_library(eigen_demo
	eigen_demo.cpp
	eigen_demo_IO.cpp
	)
