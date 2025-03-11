cmake_minimum_required(VERSION 3.10)
project(common)

add_library(COMMON_INCLUDES INTERFACE)
target_include_directories(COMMON_INCLUDES INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR}/config
        ${CMAKE_CURRENT_SOURCE_DIR}/solver_interface/include
)

file(GLOB_RECURSE COMMON_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/solver_interface/src/*.cpp)
target_sources(COMMON_INCLUDES INTERFACE ${COMMON_SOURCES})
