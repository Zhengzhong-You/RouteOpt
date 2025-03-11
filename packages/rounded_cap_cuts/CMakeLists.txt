cmake_minimum_required(VERSION 3.10)
project(rounded_cap_cuts)

add_library(ROUNDED_CAP_CUTS_INCLUDES INTERFACE)
target_include_directories(ROUNDED_CAP_CUTS_INCLUDES INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR}/chg_rc_getter/include
        ${CMAKE_CURRENT_SOURCE_DIR}/coefficient_getter/include
        ${CMAKE_CURRENT_SOURCE_DIR}/common/include
        ${CMAKE_CURRENT_SOURCE_DIR}/separation/include
)

file(GLOB_RECURSE ROUNDED_CAP_CUTS_SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/chg_rc_getter/src/*.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/coefficient_getter/src/*.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/common/src/*.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/separation/src/*.cpp
)

target_sources(ROUNDED_CAP_CUTS_INCLUDES INTERFACE ${ROUNDED_CAP_CUTS_SOURCES})