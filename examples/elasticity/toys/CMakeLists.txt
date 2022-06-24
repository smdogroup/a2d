# Add executables
add_executable(
  toy_profiler
  toy_profiler.cpp
)


# Set include dirs for the test executables
target_include_directories(toy_profiler PRIVATE ${CMAKE_SOURCE_DIR})