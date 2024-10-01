# Add targets
add_executable(test_ad_expressions test_ad_expressions.cpp)

# include A2D and test headers
target_include_directories(test_ad_expressions PRIVATE
    ${PROJECT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/tests)

# So that ctest could recognize this test
add_test(NAME test_ad_expressions COMMAND test_ad_expressions)
