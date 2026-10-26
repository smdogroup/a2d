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
