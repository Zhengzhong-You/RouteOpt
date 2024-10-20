cmake_minimum_required(VERSION 3.16)
project(CVRP)
set(CMAKE_CXX_STANDARD 17)

find_package(OpenMP REQUIRED)

option(CXX "enable C++ compilation" ON)
if(CXX)
    enable_language(CXX)
endif()

set(Boost_ROOT "/home/yzz/boost_1_83_0")
set(Boost_LIBRARY_DIR "/home/yzz/boost_1_83_0/stage/lib")

set(ENV{XGB_HOME} ./../../../Dependency/xgboost)
set(ENV{CVRPSEP_HOME} ./../../../Dependency/cvrpsep)
set(ENV{EIGEN_HOME} ./../../../Dependency/eigen-3.4.0)
#set(ENV{DELUXING_HOME} ./../../../Dependency/deLuxing)
set(ENV{HGS_HOME} ./../../../Dependency/hgs)


if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")


set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/package)
find_package(Boost REQUIRED)
find_package(GUROBI REQUIRED)
find_package(XGB REQUIRED)
find_package(CVRPSEP REQUIRED)
#find_package(DELUXING REQUIRED)
find_package(HGS REQUIRED)

include_directories(
        ${Boost_INCLUDE_DIRS}
        ${GUROBI_INCLUDE_DIRS}
        ${XGB_INCLUDE_DIRS}
        ${CVRPSEP_INCLUDE_DIRS}
#        ${DELUXING_INCLUDE_DIRS}
        $ENV{EIGEN_HOME}
        ${HGS_INCLUDE_DIRS}
)

# Automatically add all subdirectories containing header files
file(GLOB_RECURSE HEADER_LIST LIST_DIRECTORIES true "${PROJECT_SOURCE_DIR}/include/*.hpp")
foreach (HEADER ${HEADER_LIST})
    get_filename_component(HEADER_DIR ${HEADER} DIRECTORY)
    include_directories(${HEADER_DIR})
endforeach ()

file(GLOB_RECURSE SOURCES "${PROJECT_SOURCE_DIR}/src/*.cpp")


add_executable(${PROJECT_NAME} ${SOURCES})


target_link_libraries(${CMAKE_PROJECT_NAME}
        ${GUROBI_LIBRARIES}
        ${Boost_MATH_LIBRARY}
      #  ${DELUXING_LIBRARIES}
	${XGB_LIBRARIES}
       ${CVRPSEP_LIBRARIES}
${HGS_LIBRARIES}
stdc++fs
OpenMP::OpenMP_CXX
        )

if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})
    include(FeatureSummary)
    feature_summary(WHAT ALL)
endif()
