add_executable(hello hello_kokkos.cpp)

# This is the only line required to set up all compiler/linker flags
target_link_libraries(hello kokkos)
