# Specify cmake version needed
cmake_minimum_required(VERSION 3.22)

# Set project name
project(A2D VERSION 0.1)

# Set options
option(BUILD_TESTS "Compile the unit test executables" OFF)
option(BUILD_EXAMPLES "Compile the toy executables" ON)
option(BUILD_EXTENSION "Compile the pybind11 extension" ON)
option(DEBUG "a debug build without compiler optimization" OFF)

# Specify c++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Set compiler flags
if (DEBUG)
  set(CMAKE_CXX_FLAGS "-g")
else()
  set(CMAKE_CXX_FLAGS "-g -O3")
endif()

# Require Lapack
find_package(LAPACK REQUIRED)

# Require MPI
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

# Require OpenMP
find_package(OpenMP REQUIRED)

# Build pybind11 extension
if (BUILD_EXTENSION)
  execute_process(COMMAND python3 -m pybind11 --cmakedir OUTPUT_VARIABLE pybind11_DIR OUTPUT_STRIP_TRAILING_WHITESPACE)
  find_package(pybind11 REQUIRED)
  add_subdirectory(pybind11)
endif()

# Build unit tests
if (BUILD_TESTS)
  # Download google test (git tag needs to be updated frequently)
  include(FetchContent)
  FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG        e2239ee6043f73722e7aa812a459f54a28552929 # release-1.10.0
  )
  FetchContent_MakeAvailable(googletest)

  # Declare tests
  add_subdirectory(tests)
endif()

if (BUILD_EXAMPLES)
  add_subdirectory(examples/cpp)
endif()
