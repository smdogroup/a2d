# include A2D and test headers
include_directories(${A2D_ROOT_DIR}/include)
include_directories(${A2D_ROOT_DIR}/tests)

# Add targets
add_executable(test_integrand test_integrand.cpp)

# Link to kokkos
target_link_libraries(test_integrand Kokkos::kokkos)

# So that ctest could recognize this test
add_test(NAME test_integrand COMMAND test_integrand)
