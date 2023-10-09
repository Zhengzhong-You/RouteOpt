cmake_minimum_required(VERSION 3.16)
project(VRPTW)
set(CMAKE_CXX_STANDARD 17)

find_package(OpenMP REQUIRED)

option(CXX "enable C++ compilation" ON)
if(CXX)
    enable_language(CXX)
endif()

set(ENV{XGB_HOME} ./../../../Dependency/xgboost)
set(ENV{CVRPSEP_HOME} ./../../../Dependency/cvrpsep)
set(ENV{EIGEN_HOME} ./../../../Dependency/eigen-3.4.0)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")


set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/package)
find_package(GUROBI REQUIRED)
find_package(XGB REQUIRED)
find_package(CVRPSEP REQUIRED)

include_directories(${GUROBI_INCLUDE_DIRS})
include_directories(${XGB_INCLUDE_DIRS})
include_directories(${CVRPSEP_INCLUDE_DIRS})
include_directories(${PROJECT_SOURCE_DIR}/include) 
include_directories($ENV{EIGEN_HOME}) 

# list source files here
file(GLOB SOURCES "src/*.cpp")

add_executable(${PROJECT_NAME} ${SOURCES})


target_link_libraries(${CMAKE_PROJECT_NAME}
        ${GUROBI_LIBRARIES}
	${XGB_LIBRARIES}
       ${CVRPSEP_LIBRARIES}
	stdc++fs
        )

if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})
    include(FeatureSummary)
    feature_summary(WHAT ALL)
endif()
