# Add tragets
add_executable(test_matcore3d test_matcore3d.cpp)
add_executable(test_matcore2d test_matcore2d.cpp)
add_executable(test_tmp3d test_tmp3d.cpp)
add_executable(test_tmp2d test_tmp2d.cpp)

# Link to the default main from Google Test
target_link_libraries(test_matcore3d gtest_main)
target_link_libraries(test_matcore2d gtest_main)
target_link_libraries(test_tmp3d gtest_main)
target_link_libraries(test_tmp2d gtest_main)