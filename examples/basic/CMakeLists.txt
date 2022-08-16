# include A2D headers
include_directories(${A2D_ROOT_DIR}/include)

# Copy over file
configure_file(${CMAKE_SOURCE_DIR}/examples/basic/tetra_3d.vtk.in tetra_3d.vtk COPYONLY)
configure_file(${CMAKE_SOURCE_DIR}/examples/basic/tetra_3d_refine.vtk.in tetra_3d_refine.vtk COPYONLY)

# Add targets
add_executable(elasticity elasticity.cpp)
add_executable(test_grad_hess test_grad_hess.cpp)
add_executable(test_numeric test_numeric.cpp)
add_executable(toy_elasticity toy_elasticity.cpp)
add_executable(toy_helmholtz toy_helmholtz.cpp)
add_executable(toy_multiarray toy_multiarray.cpp)
add_executable(toy_profiler toy_profiler.cpp)
add_executable(toy_readvtk toy_readvtk.cpp)
add_executable(toy_tetramesh toy_tetramesh.cpp)

# Link to kokkos if used, note that linking to kokkos must happen before
# liking to OpenMP::OpenMP, otherwise it might cause compile error
if (A2D_USE_KOKKOS)
  target_link_libraries(elasticity     kokkos)
  target_link_libraries(test_grad_hess kokkos)
  target_link_libraries(test_numeric   kokkos)
  target_link_libraries(toy_elasticity kokkos)
  target_link_libraries(toy_helmholtz  kokkos)
  target_link_libraries(toy_multiarray kokkos)
  target_link_libraries(toy_profiler   kokkos)
  target_link_libraries(toy_readvtk    kokkos)
  target_link_libraries(toy_tetramesh  kokkos)
endif()

# Link libraries
target_link_libraries(elasticity     OpenMP::OpenMP_CXX LAPACK::LAPACK)
target_link_libraries(test_grad_hess OpenMP::OpenMP_CXX LAPACK::LAPACK)
target_link_libraries(test_numeric   OpenMP::OpenMP_CXX LAPACK::LAPACK)
target_link_libraries(toy_elasticity OpenMP::OpenMP_CXX LAPACK::LAPACK)
target_link_libraries(toy_helmholtz  OpenMP::OpenMP_CXX LAPACK::LAPACK)
target_link_libraries(toy_multiarray OpenMP::OpenMP_CXX LAPACK::LAPACK)
target_link_libraries(toy_profiler   OpenMP::OpenMP_CXX LAPACK::LAPACK)
target_link_libraries(toy_readvtk    OpenMP::OpenMP_CXX LAPACK::LAPACK)
target_link_libraries(toy_tetramesh  OpenMP::OpenMP_CXX LAPACK::LAPACK)

add_subdirectory(rigid_body)
