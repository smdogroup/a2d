# include A2D headers
include_directories(${A2D_ROOT_DIR}/include)

# Add targets
add_executable(parallel_element parallel_element.cpp)

# Link to kokkos, note that linking to kokkos must happen before
# liking to OpenMP::OpenMP, otherwise it might cause compile error
target_link_libraries(parallel_element Kokkos::kokkos)

# Link libraries
target_link_libraries(parallel_element OpenMP::OpenMP_CXX LAPACK::LAPACK)

# If using gcc and version < 9, need to explicitly link to filesystem
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9)
        message("Using GCC ${CMAKE_CXX_COMPILER_VERSION} < 9.0.0, explicitly link to stdc++fs")
        target_link_libraries(parallel_element stdc++fs)
    endif()
endif()
