# include A2D headers
include_directories(${A2D_ROOT_DIR}/include)

# Add targets
add_executable(elasticity elasticity.cpp)
add_executable(test_grad_hess test_grad_hess.cpp)
add_executable(test_numeric test_numeric.cpp)
add_executable(toy_elasticity toy_elasticity.cpp)
add_executable(toy_helmholtz toy_helmholtz.cpp)
add_executable(toy_multiarray toy_multiarray.cpp)
add_executable(toy_profiler toy_profiler.cpp)

# Link libraries
target_link_libraries(elasticity     OpenMP::OpenMP_CXX LAPACK::LAPACK)
target_link_libraries(test_grad_hess OpenMP::OpenMP_CXX LAPACK::LAPACK)
target_link_libraries(test_numeric   OpenMP::OpenMP_CXX LAPACK::LAPACK)
target_link_libraries(toy_elasticity OpenMP::OpenMP_CXX LAPACK::LAPACK)
target_link_libraries(toy_helmholtz  OpenMP::OpenMP_CXX LAPACK::LAPACK)
target_link_libraries(toy_multiarray OpenMP::OpenMP_CXX LAPACK::LAPACK)
target_link_libraries(toy_profiler   OpenMP::OpenMP_CXX LAPACK::LAPACK)

# Link to kokkos if used
if (A2D_USE_KOKKOS)
  target_link_libraries(elasticity     kokkos)
  target_link_libraries(test_grad_hess kokkos)
  target_link_libraries(test_numeric   kokkos)
  target_link_libraries(toy_elasticity kokkos)
  target_link_libraries(toy_helmholtz  kokkos)
  target_link_libraries(toy_multiarray kokkos)
  target_link_libraries(toy_profiler   kokkos)

  # if (Kokkos_ENABLE_CUDA)
  #   target_link_libraries(elasticity     a2dcuda)
  #   target_link_libraries(test_grad_hess a2dcuda)
  #   target_link_libraries(test_numeric   a2dcuda)
  #   target_link_libraries(toy_elasticity a2dcuda)
  #   target_link_libraries(toy_helmholtz  a2dcuda)
  #   target_link_libraries(toy_multiarray a2dcuda)
  #   target_link_libraries(toy_profiler   a2dcuda)
  # endif()

endif()