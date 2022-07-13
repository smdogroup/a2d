# Add pybind targets
pybind11_add_module(pya2d pya2d.cpp)
pybind11_add_module(pya2d_2d pya2d_2d.cpp)

# Link libraries
target_link_libraries(pya2d PUBLIC OpenMP::OpenMP_CXX LAPACK::LAPACK)
target_link_libraries(pya2d_2d PUBLIC OpenMP::OpenMP_CXX LAPACK::LAPACK)

# Copy binaries to ./a2d
install(TARGETS pya2d pya2d_2d DESTINATION ${CMAKE_SOURCE_DIR}/a2d)