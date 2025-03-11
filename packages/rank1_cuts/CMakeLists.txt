cmake_minimum_required(VERSION 3.10)
project(rank1_cuts)

add_library(RANK1_CUTS_INCLUDES INTERFACE
        coefficient_getter/src/recover_r1c_in_enum.cpp)
target_include_directories(RANK1_CUTS_INCLUDES INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR}/chg_rc_getter/include
        ${CMAKE_CURRENT_SOURCE_DIR}/coefficient_getter/include
        ${CMAKE_CURRENT_SOURCE_DIR}/common/include
        ${CMAKE_CURRENT_SOURCE_DIR}/separation/include
)

file(GLOB_RECURSE RANK1_CUTS_SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/chg_rc_getter/src/*.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/coefficient_getter/src/*.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/common/src/*.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/separation/src/*.cpp
)

target_sources(RANK1_CUTS_INCLUDES INTERFACE ${RANK1_CUTS_SOURCES})