add_executable(view view.cpp)

# This is the only line required to set up all compiler/linker flags
target_link_libraries(view kokkos)
