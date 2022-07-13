add_executable(toy_amgx toy_amgx.cpp)
configure_file(${CMAKE_SOURCE_DIR}/examples/amgx/basic/config.json config.json COPYONLY)
configure_file(${CMAKE_SOURCE_DIR}/examples/amgx/basic/matrix.mtx matrix.mtx COPYONLY)

target_link_libraries(toy_amgx PUBLIC OpenMP::OpenMP_CXX LAPACK::LAPACK)
