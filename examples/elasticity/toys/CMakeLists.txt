# Add executables
add_executable(
  toy_profiler
  toy_profiler.cpp
)

add_executable(
  toy_elasticity
  toy_elasticity.cpp
)
target_link_libraries(toy_elasticity PUBLIC OpenMP::OpenMP_CXX MPI::MPI_CXX)

add_executable(
  toy_helmholtz
  toy_helmholtz.cpp
)
target_link_libraries(toy_helmholtz PUBLIC OpenMP::OpenMP_CXX MPI::MPI_CXX)

add_executable(
  toy_multiarray
  toy_multiarray.cpp
)


# Set include dirs for the test executables
target_include_directories(toy_profiler PRIVATE ${CMAKE_SOURCE_DIR})
target_include_directories(toy_elasticity PRIVATE ${CMAKE_SOURCE_DIR})
target_include_directories(toy_helmholtz PRIVATE ${CMAKE_SOURCE_DIR})
target_include_directories(toy_multiarray PRIVATE ${CMAKE_SOURCE_DIR})