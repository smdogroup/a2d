# add_executable(linear_elastic_amgx linear_elastic_amgx.cpp)
# configure_file(${CMAKE_SOURCE_DIR}/examples/amgx/basic/config.json config.json COPYONLY)

# add_executable(elasticity elasticity.cpp)
# target_link_libraries(elasticity 
#   PUBLIC 
#     OpenMP::OpenMP_CXX 
#     LAPACK::LAPACK
#     MPI::MPI_CXX
#     CUDA::cudart_static
#     CUDA::cublas
#     CUDA::cusparse)

# add_executable(amgx amgx.cpp)
# target_link_libraries(amgx 
#   PUBLIC 
#     OpenMP::OpenMP_CXX
#     LAPACK::LAPACK
#     CUDA::cudart_static
#     CUDA::cublas
#     CUDA::cusparse)

# add_executable(amgx_mpi amgx_mpi.cpp)
# target_link_libraries(amgx_mpi 
#   PUBLIC 
#     OpenMP::OpenMP_CXX 
#     LAPACK::LAPACK
#     MPI::MPI_CXX
#     CUDA::cudart_static
#     CUDA::cublas
#     CUDA::cusparse)

add_executable(amgx_part amgx_part.cpp)
target_link_libraries(amgx_part 
  PUBLIC 
    OpenMP::OpenMP_CXX
    MPI::MPI_CXX
    LAPACK::LAPACK
    CUDA::cudart_static
    CUDA::cublas
    CUDA::cusparse)

# add_executable(amgx_mpi_part amgx_mpi_part.cpp)
# target_link_libraries(amgx_mpi_part 
#   PUBLIC 
#     OpenMP::OpenMP_CXX 
#     LAPACK::LAPACK
#     MPI::MPI_CXX
#     CUDA::cudart_static
#     CUDA::cublas
#     CUDA::cusparse)