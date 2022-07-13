add_executable(linear_elastic_amgx linear_elastic_amgx.cpp)
configure_file(${CMAKE_SOURCE_DIR}/examples/amgx/basic/config.json config.json COPYONLY)

target_link_libraries(linear_elastic_amgx PUBLIC OpenMP::OpenMP_CXX LAPACK::LAPACK)
