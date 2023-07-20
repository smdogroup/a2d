# Specify cmake version needed
cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

# Set project name, language and c++ std
project(A2D LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Make sure a build type is set explicitly
if(NOT CMAKE_BUILD_TYPE)
  message(FATAL_ERROR "CMAKE_BUILD_TYPE must be set explicitly: cmake ... -DCMAKE_BUILD_TYPE={Release,Debug}")
endif()

IF(CMAKE_BUILD_TYPE MATCHES Debug)
  message("This is a Debug build.")
ELSEIF(CMAKE_BUILD_TYPE MATCHES Release)
  message("This is a Release build.")
ELSE()
  message(FATAL_ERROR "Unknown build type CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}, should be Release or Debug")
ENDIF()

# Set paths
set(A2D_ROOT_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(KOKKOS_INSTALL_DIR ${A2D_ROOT_DIR}/installs/kokkos)
set(KOKKOS_KERNELS_INSTALL_DIR ${A2D_ROOT_DIR}/installs/kokkos-kernels)
set(Kokkos_ROOT ${KOKKOS_INSTALL_DIR}) # for a2d cmake to find kokkos cmake
set(A2D_METIS_DIR ${A2D_ROOT_DIR}/installs/metis)

# Set options
option(A2D_BUILD_EXAMPLES "Compile the a2d examples" ON)
option(A2D_BUILD_EXAMPLES_AMGX "build amgx examples that requires AMGX and CUDA" OFF)
option(A2D_BUILD_EXTENSION "Compile the pybind11 extension" OFF)
option(A2D_BUILD_UNIT_TESTS "Compile the unit test executables" OFF)

# Set warning flags
add_compile_options(
  -Wall -Wextra -Wno-unused-variable -Wno-unused-parameter -Wno-sign-compare
)

# Require Lapack and OpenMP
find_package(LAPACK REQUIRED)
find_package(OpenMP REQUIRED)
find_package(Kokkos REQUIRED)

# Build amgx examples
if(A2D_BUILD_EXAMPLES_AMGX)
  # If use AMGX, the following directories should be valid
  set(AMGX_INCLUDE_DIR ${A2D_ROOT_DIR}/extern/amgx/include)
  set(AMGX_LIB_DIR ${A2D_ROOT_DIR}/extern/amgx/lib/)

  if(NOT EXISTS ${AMGX_INCLUDE_DIR})
    message(WARNING "${AMGX_INCLUDE_DIR} does not exist, can't build with AMGX")
  elseif(NOT EXISTS ${AMGX_LIB_DIR})
    message(WARNING "${AMGX_LIB_DIR} does not exist, please build AMGX first")
  else()
    add_subdirectory(${A2D_ROOT_DIR}/examples/amgx)
  endif()
endif()

# Build pybind11 extension
if(A2D_BUILD_EXTENSION)
  execute_process(COMMAND python3 -m pybind11 --cmakedir OUTPUT_VARIABLE
    pybind11_DIR OUTPUT_STRIP_TRAILING_WHITESPACE)
  find_package(pybind11 REQUIRED)
  add_subdirectory(${A2D_ROOT_DIR}/pybind11)
endif()

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
