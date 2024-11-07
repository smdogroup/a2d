# Add targets
add_executable(test_a2dgemmcore test_a2dgemmcore.cpp)

# include A2D and test headers
target_include_directories(test_a2dgemmcore PRIVATE
    ${PROJECT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/tests)

# For tests implmented using gtest, link them to gtest
target_link_libraries(test_a2dgemmcore PRIVATE gtest_main)

include(GoogleTest)
gtest_discover_tests(test_a2dgemmcore)
