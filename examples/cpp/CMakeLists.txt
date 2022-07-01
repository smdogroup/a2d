# Add targets
add_executable(elasticity elasticity.cpp)
add_executable(test_grad_hess test_grad_hess.cpp)
add_executable(test_numeric test_numeric.cpp)
add_executable(toy_elasticity toy_elasticity.cpp)
add_executable(toy_helmholtz toy_helmholtz.cpp)
add_executable(toy_multiarray toy_multiarray.cpp)
add_executable(toy_profiler toy_profiler.cpp)

# Set include path
include_directories(${CMAKE_SOURCE_DIR}/include)

# Link libraries
target_link_libraries(elasticity     PUBLIC OpenMP::OpenMP_CXX MPI::MPI_CXX LAPACK::LAPACK)
target_link_libraries(test_grad_hess PUBLIC OpenMP::OpenMP_CXX MPI::MPI_CXX LAPACK::LAPACK)
target_link_libraries(test_numeric   PUBLIC OpenMP::OpenMP_CXX MPI::MPI_CXX LAPACK::LAPACK)
target_link_libraries(toy_elasticity PUBLIC OpenMP::OpenMP_CXX MPI::MPI_CXX LAPACK::LAPACK)
target_link_libraries(toy_helmholtz  PUBLIC OpenMP::OpenMP_CXX MPI::MPI_CXX LAPACK::LAPACK)
target_link_libraries(toy_multiarray PUBLIC OpenMP::OpenMP_CXX MPI::MPI_CXX LAPACK::LAPACK)
target_link_libraries(toy_profiler   PUBLIC OpenMP::OpenMP_CXX MPI::MPI_CXX LAPACK::LAPACK)
