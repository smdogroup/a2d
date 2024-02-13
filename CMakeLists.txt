cmake_minimum_required(VERSION 3.12)

project("A2D"
  VERSION 1.0.0
  DESCRIPTION "A light-weight header only library for 1st and 2nd order automatic differentiation"
  HOMEPAGE_URL "https://github.com/smdogroup/a2d"
  LANGUAGES CXX)

# Change the default value of CMAKE_INSTALL_PREFIX
set(A2D_INSTALL_DIR $ENV{HOME}/installs/a2d)

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX ${A2D_INSTALL_DIR} CACHE PATH "The default\
  path to install A2D to" FORCE)
  message("The directory to install A2D is set to default:\
   ${A2D_INSTALL_DIR}")
endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

option(A2D_BUILD_TESTS "Build unit tests" OFF)
option(A2D_INSTALL_LIBRARY "Enable installation" ${PROJECT_IS_TOP_LEVEL})

add_library(${PROJECT_NAME} INTERFACE)
target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_17)
target_include_directories(
  ${PROJECT_NAME}
  INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

# Set warning flags
# TODO: specify warning flags for other compilers
if(CMAKE_CXX_COMPILER_ID MATCHES "AppleClang|GNU")
  target_compile_options(${PROJECT_NAME} INTERFACE
    -Wall -Wextra -Wno-unused-variable -Wno-unused-parameter -Wno-sign-compare
  )
endif()

if(A2D_INSTALL_LIBRARY)
  install(
    TARGETS ${PROJECT_NAME}
    EXPORT ${PROJECT_NAME}_Targets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )

  write_basic_package_version_file("${PROJECT_NAME}ConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion)

  configure_package_config_file(
    "${PROJECT_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.cmake.in"
    "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
    INSTALL_DESTINATION
    ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/cmake)

  # install cmake files
  install(
    EXPORT ${PROJECT_NAME}_Targets
    FILE "${PROJECT_NAME}Targets.cmake"
    NAMESPACE ${PROJECT_NAME}::
    DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/cmake
  )
  install(
    FILES "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
    "${PROJECT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/cmake
  )

  # Install headers
  install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/ DESTINATION include)
endif()

if(A2D_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()
