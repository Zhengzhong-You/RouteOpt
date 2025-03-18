cmake_minimum_required(VERSION 3.10)
project(deluxing)

add_library(DELUXING_INCLUDES INTERFACE)
target_include_directories(DELUXING_INCLUDES INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

file(GLOB_RECURSE DELUXING_SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp
)

target_sources(DELUXING_INCLUDES INTERFACE ${DELUXING_SOURCES})