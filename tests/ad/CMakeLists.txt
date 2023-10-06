# include A2D and test headers
include_directories(${A2D_ROOT_DIR}/include)
include_directories(${A2D_ROOT_DIR}/tests)

# Add targets
add_executable(test_ad_expressions test_ad_expressions.cpp)

# Link to kokkos
target_link_libraries(test_ad_expressions Kokkos::kokkos LAPACK::LAPACK)

# So that ctest could recognize this test
add_test(NAME test_ad_expressions COMMAND test_ad_expressions)
