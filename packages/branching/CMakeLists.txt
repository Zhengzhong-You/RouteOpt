cmake_minimum_required(VERSION 3.10)
project(branching)

add_library(BRANCHING_INCLUDES INTERFACE)
target_include_directories(BRANCHING_INCLUDES INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR}/bbt/include
        ${CMAKE_CURRENT_SOURCE_DIR}/bkf/include
        ${CMAKE_CURRENT_SOURCE_DIR}/common/include
        ${CMAKE_CURRENT_SOURCE_DIR}/candidate_selector/include
)

file(GLOB_RECURSE BRANCHING_SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/bbt/src/*.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/bkf/src/*.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/common/src/*.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/candidate_selector/src/*.cpp
)

target_sources(BRANCHING_INCLUDES INTERFACE ${BRANCHING_SOURCES})