# include A2D and test headers
include_directories(${A2D_ROOT_DIR}/include)
include_directories(${A2D_ROOT_DIR}/tests)

# Add targets
add_executable(test_bsr_to_csr_csc test_bsr_to_csr_csc.cpp)

# Link to kokkos
target_link_libraries(test_bsr_to_csr_csc Kokkos::kokkos)

# Link to the default main from Google Test
target_link_libraries(test_bsr_to_csr_csc gtest_main)

# Make tests auto-testable with CMake ctest
include(GoogleTest)
gtest_discover_tests(test_bsr_to_csr_csc)
