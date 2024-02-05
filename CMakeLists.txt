# Specify cmake version needed
cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

# Set project name, language and c++ std
project(A2D LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Make sure a build type is set explicitly
if(NOT CMAKE_BUILD_TYPE)
  message("CMAKE_BUILD_TYPE not specified, set it to Release")
  set(CMAKE_BUILD_TYPE "Release")
endif()

IF(CMAKE_BUILD_TYPE MATCHES Debug)
  message("This is a Debug build.")
ELSEIF(CMAKE_BUILD_TYPE MATCHES Release)
  message("This is a Release build.")
ELSE()
  message(FATAL_ERROR "Unknown build type CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}, should be Release or Debug")
ENDIF()

# save A2D root directory
set(A2D_ROOT_DIR ${CMAKE_CURRENT_SOURCE_DIR})

# Set options
option(A2D_BUILD_EXAMPLES "Compile the a2d examples" ON)
option(A2D_BUILD_UNIT_TESTS "Compile the unit test executables" ON)

# Set warning flags
add_compile_options(
  -Wall -Wextra -Wno-unused-variable -Wno-unused-parameter -Wno-sign-compare
)

# Build unit tests
if(A2D_BUILD_UNIT_TESTS)
  # Download google test (git tag needs to be updated frequently)
  include(FetchContent)
  FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG e2239ee6043f73722e7aa812a459f54a28552929 # release-1.10.0
  )
  FetchContent_MakeAvailable(googletest)
  add_subdirectory(${A2D_ROOT_DIR}/tests)
endif()

# Build a2d examples
if(A2D_BUILD_EXAMPLES)
  add_subdirectory(${A2D_ROOT_DIR}/examples)
endif()
