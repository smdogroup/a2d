# Include A2D and AMGX headers
include_directories(${AMGX_INCLUDE_DIR})
include_directories(${A2D_ROOT_DIR}/include)
include_directories(${A2D_ROOT_DIR}/include/cuda)
include_directories(${A2D_ROOT_DIR}/include/sparse)

link_directories(${AMGX_LIB_DIR})
link_libraries(${AMGX_LIB_DIR}/libamgxsh.so)

configure_file(${CMAKE_SOURCE_DIR}/examples/amgx/config.json config.json COPYONLY)
configure_file(${CMAKE_SOURCE_DIR}/examples/amgx/matrix.mtx matrix.mtx COPYONLY)

add_executable(amgx amgx.cpp)
target_link_libraries(amgx 
  PUBLIC 
    OpenMP::OpenMP_CXX
    LAPACK::LAPACK
    CUDA::cudart_static)

add_executable(toy_amgx toy_amgx.cpp)
target_link_libraries(toy_amgx 
  PRIVATE 
    OpenMP::OpenMP_CXX 
    LAPACK::LAPACK)