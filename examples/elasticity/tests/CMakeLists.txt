# Add executables
add_executable(
  test_matcore
  test_matcore.cpp
)

add_executable(
  test_matcore2d
  test_matcore2d.cpp
)

add_executable(
  test_tmp
  test_tmp.cpp
)

add_executable(
  test_tmp2d
  test_tmp2d.cpp
)

# Set include dirs for the test executables
target_include_directories(test_matcore PRIVATE ${CMAKE_SOURCE_DIR})
target_include_directories(test_matcore2d PRIVATE ${CMAKE_SOURCE_DIR})
target_include_directories(test_tmp PRIVATE ${CMAKE_SOURCE_DIR})
target_include_directories(test_tmp2d PRIVATE ${CMAKE_SOURCE_DIR})