# FindCVRPSEP.cmake

set(CVRPSEP_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/cvrpsep")

if (EXISTS "${CVRPSEP_ROOT}")
    set(CVRPSEP_FOUND TRUE)
    message(STATUS "Found CVRPSEP root at: ${CVRPSEP_ROOT}")
else ()
    set(CVRPSEP_FOUND FALSE)
    message(FATAL_ERROR "CVRPSEP not found at ${CVRPSEP_ROOT}")
endif ()

# Directly include all headers under cvrpsep/
set(CVRPSEP_INCLUDE_DIRS "${CVRPSEP_ROOT}")

# Locate the library in legacy obj/ (Makefile) or cmake build output.
find_library(CVRPSEP_LIBRARY
        NAMES cvrpsep libcvrpsep
        PATHS
        "${CVRPSEP_ROOT}/obj"
        "${CVRPSEP_ROOT}/build"
        "${CVRPSEP_ROOT}/build/Release"
        "${CVRPSEP_ROOT}/build/RelWithDebInfo"
        NO_DEFAULT_PATH
)

if (CVRPSEP_LIBRARY)
    message(STATUS "Found CVRPSEP library: ${CVRPSEP_LIBRARY}")
else ()
    message(FATAL_ERROR "CVRPSEP library not found in ${CVRPSEP_ROOT}/obj")
endif ()

set(CVRPSEP_LIBRARIES "${CVRPSEP_LIBRARY}")

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CVRPSEP DEFAULT_MSG CVRPSEP_LIBRARY CVRPSEP_INCLUDE_DIRS)

mark_as_advanced(CVRPSEP_INCLUDE_DIRS CVRPSEP_LIBRARY)
