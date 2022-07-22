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
target_link_libraries(elasticity     PRIVATE OpenMP::OpenMP_CXX LAPACK::LAPACK)
target_link_libraries(test_grad_hess PRIVATE OpenMP::OpenMP_CXX LAPACK::LAPACK)
target_link_libraries(test_numeric   PRIVATE OpenMP::OpenMP_CXX LAPACK::LAPACK)
target_link_libraries(toy_elasticity PRIVATE OpenMP::OpenMP_CXX LAPACK::LAPACK)
target_link_libraries(toy_helmholtz  PRIVATE OpenMP::OpenMP_CXX LAPACK::LAPACK)
target_link_libraries(toy_multiarray PRIVATE OpenMP::OpenMP_CXX LAPACK::LAPACK)
target_link_libraries(toy_profiler   PRIVATE OpenMP::OpenMP_CXX LAPACK::LAPACK)