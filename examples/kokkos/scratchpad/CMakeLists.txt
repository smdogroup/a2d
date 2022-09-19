add_executable(scratchpad scratchpad.cpp)

# This is the only line required to set up all compiler/linker flags
target_link_libraries(scratchpad kokkos)
