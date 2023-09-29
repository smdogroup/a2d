# include A2D headers
include_directories(${A2D_ROOT_DIR}/include)

# Set up paropt library and dependencies
find_package(MPI REQUIRED)
link_directories($ENV{HOME}/git/paropt/lib)
include_directories(${MPI_INCLUDE_PATH})
include_directories($ENV{HOME}/git/paropt/src)

# include metis
include_directories(${A2D_METIS_DIR}/include)
link_directories(${A2D_METIS_DIR}/lib)

# Add targets
add_executable(2d_analysis 2d_analysis.cpp)
add_executable(2d_topology 2d_topology.cpp)

# Link to kokkos, note that linking to kokkos must happen before
# liking to OpenMP::OpenMP, otherwise it might cause compile error
target_link_libraries(2d_analysis Kokkos::kokkos)
target_link_libraries(2d_topology Kokkos::kokkos)

# Link libraries
target_link_libraries(2d_analysis OpenMP::OpenMP_CXX LAPACK::LAPACK metis)
target_link_libraries(2d_topology OpenMP::OpenMP_CXX LAPACK::LAPACK metis)

target_link_libraries(2d_topology paropt MPI::MPI_CXX)
target_link_libraries(2d_analysis paropt MPI::MPI_CXX)

# If using gcc and version < 9, need to explicitly link to filesystem
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9)
        message("Using GCC ${CMAKE_CXX_COMPILER_VERSION} < 9.0.0, explicitly link to stdc++fs")
        target_link_libraries(2d_analysis stdc++fs)
        target_link_libraries(2d_topology stdc++fs)
    endif()
endif()
