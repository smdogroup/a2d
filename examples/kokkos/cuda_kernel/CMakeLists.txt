add_executable(cuda_kernel cuda_kernel.cpp)

# This is the only line required to set up all compiler/linker flags
target_link_libraries(cuda_kernel kokkos)
