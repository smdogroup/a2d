# include A2D headers
include_directories(${A2D_ROOT_DIR}/include)

# Copy over file
configure_file(${CMAKE_SOURCE_DIR}/examples/test_verify/tetra_3d.vtk.in tetra_3d.vtk COPYONLY)

# Add targets
add_executable(verify_element verify_element.cpp)
add_executable(verify_L2element verify_L2element.cpp)
add_executable(test_profiler test_profiler.cpp)
add_executable(test_multiphysics test_multiphysics.cpp)
add_executable(test_readvtk test_readvtk.cpp)

# Link to kokkos, note that linking to kokkos must happen before
# liking to OpenMP::OpenMP, otherwise it might cause compile error
target_link_libraries(verify_element Kokkos::kokkos)
target_link_libraries(verify_L2element Kokkos::kokkos)
target_link_libraries(test_profiler Kokkos::kokkos)
target_link_libraries(test_readvtk Kokkos::kokkos)
target_link_libraries(test_multiphysics Kokkos::kokkos)

# Link libraries
target_link_libraries(verify_element OpenMP::OpenMP_CXX LAPACK::LAPACK)
target_link_libraries(verify_L2element OpenMP::OpenMP_CXX LAPACK::LAPACK)
target_link_libraries(test_profiler OpenMP::OpenMP_CXX LAPACK::LAPACK)
target_link_libraries(test_readvtk OpenMP::OpenMP_CXX LAPACK::LAPACK)
target_link_libraries(test_multiphysics OpenMP::OpenMP_CXX LAPACK::LAPACK)

# If using gcc and version < 9, need to explicitly link to filesystem
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9)
        message("Using GCC ${CMAKE_CXX_COMPILER_VERSION} < 9.0.0, explicitly link to stdc++fs")
        target_link_libraries(verify_element stdc++fs)
        target_link_libraries(verify_L2element stdc++fs)
        target_link_libraries(test_profiler stdc++fs)
        target_link_libraries(test_readvtk stdc++fs)
        target_link_libraries(test_multiphysics stdc++fs)
    endif()
endif()
