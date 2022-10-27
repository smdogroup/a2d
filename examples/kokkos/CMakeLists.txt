# Include A2D headers
include_directories(${A2D_ROOT_DIR}/include)

# Which .cpp to compile
add_executable(scratchpad scratchpad.cpp)

# Link libs
target_link_libraries(scratchpad Kokkos::kokkos)