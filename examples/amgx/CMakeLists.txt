include_directories(${AMGX_INCLUDE_DIR})
link_directories(${AMGX_LIB_DIR})
link_libraries(${AMGX_LIB_DIR}/libamgxsh.so)

add_subdirectory(basic)
add_subdirectory(linear_elastic)