# include A2D and test headers
include_directories(${A2D_ROOT_DIR}/include)
include_directories(${A2D_ROOT_DIR}/tests)
include_directories(${A2D_METIS_DIR}/include)

# link to metis
link_directories(${A2D_METIS_DIR}/lib)

# Add targets
add_executable(test_analysis test_analysis.cpp)

# target_compile_options(test_analysis PRIVATE -fsanitize=address)
# target_link_options(test_analysis PRIVATE -fsanitize=address)

# Link to kokkos
target_link_libraries(test_analysis Kokkos::kokkos OpenMP::OpenMP_CXX LAPACK::LAPACK metis)

# So that ctest could recognize this test
add_test(NAME test_analysis COMMAND test_analysis)
