set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Download google test (git tag needs to be updated frequently)
include(FetchContent)
FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG 6910c9d9165801d8827d628cb72eb7ea9dd538c5 # release-1.16.0
)
FetchContent_MakeAvailable(googletest)

# Add subdirectories
add_subdirectory(ad)

# Add individual tests
add_executable(test_a2dtuple test_a2dtuple.cpp)

# include A2D and test headers
target_include_directories(test_a2dtuple PRIVATE
    ${PROJECT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/tests)

# For tests implmented using gtest, link them to gtest
target_link_libraries(test_a2dtuple PRIVATE gtest_main)

include(GoogleTest)
gtest_discover_tests(test_a2dtuple)
