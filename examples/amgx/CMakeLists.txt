# Include A2D and AMGX headers
include_directories(${AMGX_INCLUDE_DIR})
include_directories(${A2D_ROOT_DIR}/include)

link_directories(${AMGX_LIB_DIR})
link_libraries(${AMGX_LIB_DIR}/libamgxsh.so)

add_subdirectory(basic)
add_subdirectory(linear_elastic)