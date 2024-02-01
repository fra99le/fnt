# CMakeLists.txt
# fnt: Numeric Toolbox
#
# Copyright (c) 2024 Bryan Franklin. All rights reserved.
#
cmake_minimum_required(VERSION 3.0)
project (libfnt)

# remove the 'lib' prefix from modules
set(CMAKE_SHARED_MODULE_PREFIX "")

if(${CMAKE_C_COMPILER_ID} MATCHES "Clang" )
    if( CMAKE_MAJOR_VERSION GREATER_EQUAL 3 AND CMAKE_MINOR_VERSION GREATER_EQUAL 13 )
        add_link_options("SHELL:-undefined dynamic_lookup")
    else()
        set(CMAKE_MODULE_LINKER_FLAGS "-undefined dynamic_lookup")
    endif()
endif()

set(SOURCE_FILES "")
AUX_SOURCE_DIRECTORY(. SOURCE_FILES)
foreach(SRC ${SOURCE_FILES})
    get_filename_component(BASE ${SRC} NAME_WE)
    message(BASE="${BASE}")
    add_library(${BASE} MODULE ${SRC})
    set_property(TARGET ${BASE} PROPERTY C_STANDARD 99)
endforeach(SRC)
