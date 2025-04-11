cmake_minimum_required(VERSION 3.10)
project(cvrp VERSION 3.0)


set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif ()

option(ENABLE_VALGRIND_MEM_CHECK "Enable Valgrind Memory Check" OFF)

if (ENABLE_VALGRIND_MEM_CHECK)
    add_compile_definitions(VALGRIND_MEM_CHECK)
endif ()


set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../../external external_build)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../../common common_build)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../../branching branching_build)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../../rounded_cap_cuts rcc_build)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../../rank1_cuts rank1_cuts_build)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../../deluxing deluxing_build)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/lib/cmake_modules")
find_package(HGS REQUIRED)
find_package(ZLIB REQUIRED)

# Find the required packages of the external
include_directories(
        ${BOOST_INCLUDE_DIRS} # no lib
        ${EIGEN_INCLUDE_DIRS} # no lib
        ${CVRPSEP_INCLUDE_DIRS}
        ${GUROBI_INCLUDE_DIRS}
        ${XGB_INCLUDE_DIRS}
        ${HGS_INCLUDE_DIRS}
)


# Include directories
include_directories(
        ${CMAKE_CURRENT_SOURCE_DIR}/src/main/include
        ${CMAKE_CURRENT_SOURCE_DIR}/src/vrptw/include
        ${CMAKE_CURRENT_SOURCE_DIR}/src/add_column/include
        ${CMAKE_CURRENT_SOURCE_DIR}/src/cuts/include
        ${CMAKE_CURRENT_SOURCE_DIR}/src/pricing/include
        ${CMAKE_CURRENT_SOURCE_DIR}/src/l2b/include
        ${CMAKE_CURRENT_SOURCE_DIR}/src/read_data/include
        ${CMAKE_CURRENT_SOURCE_DIR}/src/heuristics/include
        ${CMAKE_CURRENT_SOURCE_DIR}/src/initial/include
        ${CMAKE_CURRENT_SOURCE_DIR}/src/node/include
        ${CMAKE_CURRENT_SOURCE_DIR}/src/two_stage/include
)


file(GLOB_RECURSE CVRP_SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/src/main/src/*.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/vrptw/src/*.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/add_column/src/*.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/cuts/src/*.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/pricing/src/*.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/l2b/src/*.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/read_data/src/*.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/heuristics/src/*.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/initial/src/*.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/node/src/*.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/two_stage/src/*.cpp
)


add_executable(${PROJECT_NAME} ${CVRP_SOURCES}
)

# set_target_properties(${PROJECT_NAME} PROPERTIES OUTPUT_NAME "cvrp_rootcutonly")


target_link_libraries(${PROJECT_NAME}
        PRIVATE
        COMMON_INCLUDES
        BRANCHING_INCLUDES
        ROUNDED_CAP_CUTS_INCLUDES
        RANK1_CUTS_INCLUDES
        DELUXING_INCLUDES
        ${GUROBI_LIBRARIES}
        ${CVRPSEP_LIBRARIES}
        ${XGB_LIBRARIES}
        ${HGS_LIBRARIES}
        ZLIB::ZLIB
)
