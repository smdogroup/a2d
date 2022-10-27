# include A2D headers
include_directories(${A2D_ROOT_DIR}/include)

# Copy over the mesh
if(EXISTS ${CMAKE_SOURCE_DIR}/extern/mesh/rigid_body.vtk)
  configure_file(${CMAKE_SOURCE_DIR}/extern/mesh/rigid_body.vtk rigid_body.vtk COPYONLY)
else()
  message(WARNING "${CMAKE_SOURCE_DIR}/extern/mesh/rigid_body.vtk does not exist")
endif()

# Add target
add_executable(rigid_body rigid_body.cpp)

# Link
target_link_libraries(rigid_body Kokkos::kokkos OpenMP::OpenMP_CXX LAPACK::LAPACK)
