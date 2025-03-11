# FindGUROBI.cmake

set(GUROBI_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/gurobi1201/linux64")

if (EXISTS "${GUROBI_ROOT}")
    set(GUROBI_FOUND TRUE)
    message(STATUS "Found GUROBI root at: ${GUROBI_ROOT}")
else ()
    set(GUROBI_FOUND FALSE)
    message(FATAL_ERROR "GUROBI not found at ${GUROBI_ROOT}")
endif ()

find_path(GUROBI_INCLUDE_DIR
        NAMES gurobi_c.h
        PATHS "${GUROBI_ROOT}/include")

find_library(GUROBI_LIBRARY
        NAMES libgurobi120.so
        PATHS "${GUROBI_ROOT}/lib"
)

set(GUROBI_INCLUDE_DIRS "${GUROBI_INCLUDE_DIR}")
set(GUROBI_LIBRARIES "${GUROBI_LIBRARY}")

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GUROBI DEFAULT_MSG GUROBI_LIBRARY GUROBI_INCLUDE_DIR)


mark_as_advanced(GUROBI_INCLUDE_DIR GUROBI_LIBRARY)