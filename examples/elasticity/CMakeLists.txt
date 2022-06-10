# Specify cmake version needed
cmake_minimum_required(VERSION 3.22)

# Set project name
project(A2D VERSION 0.1)

# Specify c++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Set compiler flags
set(CMAKE_CXX_FLAGS "-g -O3 -Wall -Wextra")

# Require MPI
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

# Require OpenMP
find_package(OpenMP REQUIRED)

# Build executables
add_executable(test_grad_hess test_grad_hess.cpp)
target_link_libraries(test_grad_hess PUBLIC OpenMP::OpenMP_CXX MPI::MPI_CXX)

add_executable(elasticity elasticity.cpp)
target_link_libraries(elasticity PUBLIC OpenMP::OpenMP_CXX MPI::MPI_CXX)

# Require pybind11 to build the extension
execute_process(COMMAND python3 -m pybind11 --cmakedir OUTPUT_VARIABLE pybind11_DIR OUTPUT_STRIP_TRAILING_WHITESPACE)
find_package(pybind11 REQUIRED)

# Build extension
pybind11_add_module(example example.cpp)

# Download google test (commit id needs to be updated frequently)
include(FetchContent)
FetchContent_Declare(
   googletest
   URL https://github.com/google/googletest/archive/b1f84bf1763b1010597bff13c79b5388eebdf205.zip
)
FetchContent_MakeAvailable(googletest)

# Declare tests
add_subdirectory(tests)

# Link to the default main from Google Test
target_link_libraries(
  test_matcore
  gtest_main
)
