include_directories(
  ${CMAKE_SOURCE_DIR}/ext
  ${CMAKE_SOURCE_DIR}/src
)

#add_subdirectory(stiffness_checker)

add_subdirectory(eigen_demo)
add_subdirectory(bindings/pyconmech)
