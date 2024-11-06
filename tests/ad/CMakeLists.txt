# Add targets
add_executable(test_ad_expressions test_ad_expressions.cpp)
add_executable(test_a2dmat test_a2dmat.cpp)
add_executable(test_a2dgemm test_a2dgemm.cpp)

# include A2D and test headers
target_include_directories(test_ad_expressions PRIVATE
    ${PROJECT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/tests)
target_include_directories(test_a2dmat PRIVATE
    ${PROJECT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/tests)
target_include_directories(test_a2dgemm PRIVATE
    ${PROJECT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/tests)

# For tests implmented using gtest, link them to gtest
target_link_libraries(test_a2dmat PRIVATE gtest_main)
target_link_libraries(test_a2dgemm PRIVATE gtest_main)

include(GoogleTest)
gtest_discover_tests(test_a2dmat)
gtest_discover_tests(test_a2dgemm)

# Add non-gtest tests manually so that ctest could recognize it's a test
add_test(NAME test_ad_expressions COMMAND test_ad_expressions)
