# Specify cmake version needed
cmake_minimum_required(VERSION 3.22)

# Set project name
project(A2D VERSION 0.1)

# Set options
option(BUILD_TESTS "Compile the unit test executables" OFF)
option(BUILD_EXAMPLES "Compile the toy executables" ON)
option(BUILD_AMGX_EXAMPLES "build amgx examples that requires AMGX and CUDA" ON)
option(BUILD_EXTENION "Compile the pybind11 extension" ON)
option(DEBUG "a debug build without compiler optimization" OFF)

# Specify c++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Set compiler flags
if (DEBUG)
  set(CMAKE_CXX_FLAGS "-g -Wall -Wno-unused-variable -Wno-unused-parameter -Wno-sign-compare")
else()
  set(CMAKE_CXX_FLAGS "-g -O3 -Wall -Wno-unused-variable -Wno-unused-parameter -Wno-sign-compare")
endif()

# Set A2D include directory
include_directories(${CMAKE_SOURCE_DIR}/include)

# Require Lapack
find_package(LAPACK REQUIRED)

# Require MPI
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

# Require OpenMP
find_package(OpenMP REQUIRED)

# Check if AMGX is available
set(AMGX_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/extern/AMGX/include)
set(AMGX_LIB_DIR ${CMAKE_SOURCE_DIR}/extern/AMGX/lib/)

if (BUILD_AMGX_EXAMPLES)
  if (NOT EXISTS "${AMGX_INCLUDE_DIR}")
    message(WARNING "${AMGX_INCLUDE_DIR} does not exist, can't build with AMGX")
  elseif(NOT EXISTS "${AMGX_LIB_DIR}")
    message(WARNING "${AMGX_LIB_DIR} does not exist, please build AMGX first")
  else()
    add_subdirectory(examples/amgx)
  endif()
endif()

# Build pybind11 extension
if (BUILD_EXTENION)
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

# Build cpp examples
if (BUILD_EXAMPLES)
  add_subdirectory(examples/cpp)
endif()
