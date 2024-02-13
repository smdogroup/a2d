cmake_minimum_required(VERSION 3.12)

project(A2D_EXAMPLES)

if(NOT DEFINED A2D_INSTALL_DIR)
    set(A2D_INSTALL_DIR $ENV{HOME}/installs/a2d)
endif()

if(NOT EXISTS ${A2D_INSTALL_DIR})
    message(FATAL_ERROR "${A2D_INSTALL_DIR} does not exist, specify the\
    path to the A2D installation by -DA2D_INSTALL_DIR=..."
    )
endif()

find_package(A2D REQUIRED PATHS ${A2D_INSTALL_DIR})

add_subdirectory(ad)
