#
# CMakeLists.txt
# fnt: Numerical Toolbox
#
# Copyright (c) 2024 Bryan Franklin. All rights reserved.
#
cmake_minimum_required(VERSION 3.0)
project(libfnt)

# The version number.
set (FNT_VERSION_MAJOR 0)
set (FNT_VERSION_MINOR 0)
set (FNT_VERSION_PATCH 1)

# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
  "${PROJECT_SOURCE_DIR}/libfnt_config.h.in"
  "${PROJECT_BINARY_DIR}/libfnt_config.h"
  )

set (CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# set library as output
AUX_SOURCE_DIRECTORY(. SOURCE_FILES)
add_library(libfnt SHARED ${SOURCE_FILES})

# add common system libraries
target_link_libraries(libfnt m)

# check for optional libraries
#find_package(IPP_CORE)
#find_package(IPP_VM)
#find_package(IPP_S)

add_subdirectory(methods)
add_subdirectory(tests)

# add a target to remove files cmake creates
add_custom_target(full-clean
                  COMMAND bash -c 'rm -frv CMakeCache.txt CMakeFiles Makefile; rm -fvf CMake*.log'
                  )
add_custom_target(fix-cmake
                  DEPENDS full-clean
                  )

# add a debug build target
add_custom_target(debug-build
                  COMMAND bash -c 'cmake -DCMAKE_BUILD_TYPE=DEBUG -DCMAKE_C_FLAGS_DEBUG="-O0" . && CFLAGS=-Wall VERBOSE=1 make'
                  )

