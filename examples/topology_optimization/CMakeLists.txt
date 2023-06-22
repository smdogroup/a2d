# include A2D headers
include_directories(${A2D_ROOT_DIR}/include)

# Set up paropt library and dependencies
find_package(MPI REQUIRED)
link_directories($ENV{HOME}/git/paropt/lib)
include_directories(${MPI_INCLUDE_PATH})
include_directories($ENV{HOME}/git/paropt/src)

# Add targets
add_executable(toy_elasticity toy_elasticity.cpp)
add_executable(toy_grad_check toy_grad_check.cpp)
add_executable(topology_elasticity topology_elasticity.cpp)
add_executable(topology_heat topology_heat.cpp)

# Link to kokkos, note that linking to kokkos must happen before
# liking to OpenMP::OpenMP, otherwise it might cause compile error
target_link_libraries(toy_elasticity Kokkos::kokkos)
target_link_libraries(toy_grad_check Kokkos::kokkos)
target_link_libraries(topology_elasticity Kokkos::kokkos)
target_link_libraries(topology_heat Kokkos::kokkos)

# Link libraries
target_link_libraries(toy_elasticity OpenMP::OpenMP_CXX LAPACK::LAPACK)
target_link_libraries(toy_grad_check OpenMP::OpenMP_CXX LAPACK::LAPACK)
target_link_libraries(topology_elasticity OpenMP::OpenMP_CXX LAPACK::LAPACK)
target_link_libraries(topology_heat OpenMP::OpenMP_CXX LAPACK::LAPACK)

# Link optimization applications to paropt and MPI
target_link_libraries(topology_elasticity paropt MPI::MPI_CXX)
target_link_libraries(topology_heat paropt MPI::MPI_CXX)

# If using gcc and version < 9, need to explicitly link to filesystem
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9)
        message("Using GCC ${CMAKE_CXX_COMPILER_VERSION} < 9.0.0, explicitly link to stdc++fs")
        target_link_libraries(toy_elasticity stdc++fs)
        target_link_libraries(toy_grad_check stdc++fs)
        target_link_libraries(topology_elasticity stdc++fs)
        target_link_libraries(topology_heat stdc++fs)
    endif()
endif()