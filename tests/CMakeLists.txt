enable_testing()

# include A2D headers
include_directories(${A2D_ROOT_DIR}/include)
include_directories(${CMAKE_SOURCE_DIR}/include)

# Add targets
add_executable(test_matcore3d test_matcore3d.cpp)
add_executable(test_matcore2d test_matcore2d.cpp)
add_executable(test_matops3d test_matops3d.cpp)
add_executable(test_matops2d test_matops2d.cpp)
add_executable(test_vecops3d test_vecops3d.cpp)
add_executable(test_vecops3d_jg test_vecops3d_jg.cpp)
add_executable(test_vecops3d_a2d test_vecops3d_a2d.cpp)
add_executable(test_matvecops3d_ad test_matvecops3d_ad.cpp)
add_executable(test_matvecops3d_a2d test_matvecops3d_a2d.cpp)

# Link to the default main from Google Test
target_link_libraries(test_matcore3d gtest_main)
target_link_libraries(test_matcore2d gtest_main)
target_link_libraries(test_matops3d gtest_main)
target_link_libraries(test_matops2d gtest_main)
target_link_libraries(test_vecops3d gtest_main)
target_link_libraries(test_vecops3d_jg gtest_main)
target_link_libraries(test_vecops3d_a2d gtest_main)
target_link_libraries(test_matvecops3d_ad gtest_main)
target_link_libraries(test_matvecops3d_a2d gtest_main)

# Make tests auto-testable with CMake ctest
include(GoogleTest)
gtest_discover_tests(test_matcore3d)
gtest_discover_tests(test_matcore2d)
gtest_discover_tests(test_matops3d)
gtest_discover_tests(test_matops2d)
gtest_discover_tests(test_vecops3d)
gtest_discover_tests(test_vecops3d_jg)
gtest_discover_tests(test_vecops3d_a2d)
gtest_discover_tests(test_matvecops3d_ad)
gtest_discover_tests(test_matvecops3d_a2d)