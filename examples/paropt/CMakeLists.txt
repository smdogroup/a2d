# include A2D headers
include_directories(${A2D_ROOT_DIR}/include)


# Include MPI and ParOpt headers
include_directories(${MPI_INCLUDE_PATH})
include_directories($ENV{HOME}/git/paropt/src)

# Load MPI and paropt libraries
find_package(MPI REQUIRED)
link_directories($ENV{HOME}/git/paropt/lib)

# Add target and link to MPI and paropt
add_executable(paropt_example paropt_example.cpp)
target_link_libraries(paropt_example PRIVATE paropt MPI::MPI_CXX)
