# Specify cmake version needed
cmake_minimum_required(VERSION 3.22 FATAL_ERROR)

# Set project name, language and c++ std
project(A2D LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set A2D root directory
set(A2D_ROOT_DIR ${CMAKE_CURRENT_SOURCE_DIR})

# Set options
option(A2D_BUILD_UNIT_TESTS "Compile the unit test executables" OFF)
option(A2D_BUILD_EXAMPLES_BASIC "Compile the a2d toy executables" ON)
option(A2D_BUILD_EXAMPLES_AMGX "build amgx examples that requires AMGX and CUDA" OFF)
option(A2D_BUILD_EXAMPLES_KOKKOS "build kokkos examples" OFF)
option(A2D_BUILD_EXTENSION "Compile the pybind11 extension" OFF)

# Set global compiler flags
add_compile_options(
  -O3 -g -Wall -Wextra -Wno-unused-variable -Wno-unused-parameter -Wno-sign-compare
  )

# Require Lapack and OpenMP
find_package(LAPACK REQUIRED)
find_package(OpenMP REQUIRED)

# Build amgx examples
if (A2D_BUILD_EXAMPLES_AMGX)
  # If use AMGX, the following directories should be valid
  set(AMGX_INCLUDE_DIR ${A2D_ROOT_DIR}/extern/amgx/include)
  set(AMGX_LIB_DIR ${A2D_ROOT_DIR}/extern/amgx/lib/)
  if (NOT EXISTS ${AMGX_INCLUDE_DIR})
    message(WARNING "${AMGX_INCLUDE_DIR} does not exist, can't build with AMGX")
  elseif(NOT EXISTS ${AMGX_LIB_DIR})
    message(WARNING "${AMGX_LIB_DIR} does not exist, please build AMGX first")
  else()
    add_subdirectory(${A2D_ROOT_DIR}/examples/amgx)
  endif()
endif()

# If use kokkos, we build kokkos in-place
set(A2D_KOKKOS_DIR ${A2D_ROOT_DIR}/extern/kokkos)
if (NOT EXISTS ${A2D_KOKKOS_DIR})
  message(FATAL_ERROR "${A2D_KOKKOS_DIR} does not exist, can't use Kokkos")
else()
  # Build Kokkos
  add_subdirectory(${A2D_KOKKOS_DIR})
  include_directories(${A2D_KOKKOS_DIR}/core/src)
  include_directories(${A2D_KOKKOS_DIR}/algorithms/src)
  include_directories(${A2D_KOKKOS_DIR}/containers/src)
  include_directories(${CMAKE_BINARY_DIR}/extern/kokkos)
  include_directories(${CMAKE_BINARY_DIR}/extern/kokkos/algorithms/src)
  include_directories(${CMAKE_BINARY_DIR}/extern/kokkos/containers/src)
endif()

# Build pybind11 extension
if (A2D_BUILD_EXTENSION)
  execute_process(COMMAND python3 -m pybind11 --cmakedir OUTPUT_VARIABLE
  pybind11_DIR OUTPUT_STRIP_TRAILING_WHITESPACE)
  find_package(pybind11 REQUIRED)
  add_subdirectory(${A2D_ROOT_DIR}/pybind11)
endif()

# Build unit tests
if (A2D_BUILD_UNIT_TESTS)
  # Download google test (git tag needs to be updated frequently)
  include(FetchContent)
  FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG        e2239ee6043f73722e7aa812a459f54a28552929 # release-1.10.0
  )
  FetchContent_MakeAvailable(googletest)
  add_subdirectory(${A2D_ROOT_DIR}/tests)
endif()

# Build pure a2d examples
if (A2D_BUILD_EXAMPLES_BASIC)
  add_subdirectory(${A2D_ROOT_DIR}/examples/basic)
endif()

# Build kokkos examples
if (A2D_BUILD_EXAMPLES_KOKKOS)
  add_subdirectory(${A2D_ROOT_DIR}/examples/kokkos)
endif()
