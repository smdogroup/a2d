# Specify cmake version needed
cmake_minimum_required(VERSION 3.22)

# Set project name
project(A2D VERSION 0.1)

# Set options
option(BUILD_TESTS "Compile the unit test executables" OFF)
option(BUILD_TOYS "Compile the toy executables" ON)
option(DEBUG "a debug build" OFF)

# Specify c++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Set compiler flags
if (DEBUG)
  set(CMAKE_CXX_FLAGS "-g")
else()
  set(CMAKE_CXX_FLAGS "-g -O3")
endif()

# Require MPI
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

# Require OpenMP
find_package(OpenMP REQUIRED)

# Build executables
add_executable(test_grad_hess test_grad_hess.cpp)
target_link_libraries(test_grad_hess PUBLIC OpenMP::OpenMP_CXX MPI::MPI_CXX)

add_executable(test_numeric test_numeric.cpp)
target_link_libraries(test_numeric PUBLIC OpenMP::OpenMP_CXX MPI::MPI_CXX)

add_executable(elasticity elasticity.cpp)
target_link_libraries(elasticity PUBLIC OpenMP::OpenMP_CXX MPI::MPI_CXX)

# Require pybind11 to build the extension
execute_process(COMMAND python3 -m pybind11 --cmakedir OUTPUT_VARIABLE pybind11_DIR OUTPUT_STRIP_TRAILING_WHITESPACE)
find_package(pybind11 REQUIRED)

# Build extension
pybind11_add_module(example example.cpp)

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

  # Link to the default main from Google Test
  target_link_libraries(
    test_matcore
    gtest_main
  )

  target_link_libraries(
    test_matcore2d
    gtest_main
  )

  target_link_libraries(
    test_tmp
    gtest_main
  )

  target_link_libraries(
    test_tmp2d
    gtest_main
  )
endif()

if (BUILD_TOYS)
  add_subdirectory(toys)
endif()