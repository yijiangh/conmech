include_directories(
	${Eigen_INCLUDE_DIRS}
	${JSON_INCLUDE_DIRS}
)

if (NOT APPLE)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fvisibility=hidden")
endif()

# TODO: add tests
#add_subdirectory(test)

#add_library(stiffness_checker
#
#	)

#target_link_libraries(stiffness_checker )