cmake_minimum_required(VERSION 3.12)
project(hgscvrp)

set(CMAKE_CXX_STANDARD 11)

set(SOURCES
    src/Split.cpp
    src/Population.cpp
    src/Params.cpp
    src/my_hgs.cpp
    src/LocalSearch.cpp
    src/InstanceCVRPLIB.cpp
    src/Individual.cpp
    src/Genetic.cpp
    src/C_Interface.cpp
    src/AlgorithmParameters.cpp
)


include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)


add_library(hgscvrp SHARED ${SOURCES})


set_target_properties(hgscvrp PROPERTIES OUTPUT_NAME "hgscvrp")

install(TARGETS hgscvrp DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/lib)
