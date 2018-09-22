
pybind11_add_module(eigen_demo_py
	eigen_demo_py.cpp)

target_link_libraries(eigen_demo_py PRIVATE eigen_demo)
