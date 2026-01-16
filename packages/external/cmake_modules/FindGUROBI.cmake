# FindGUROBI.cmake

set(GUROBI_ROOT "/Library/gurobi1300/macos_universal2")

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

set(_GUROBI_LIB_PATTERNS
        "${GUROBI_ROOT}/lib/libgurobi[0-9]*.so"
        "${GUROBI_ROOT}/lib/libgurobi[0-9]*.dylib"
        "${GUROBI_ROOT}/lib/libgurobi[0-9]*.a"
)

file(GLOB _GUROBI_LIBS ${_GUROBI_LIB_PATTERNS})
list(FILTER _GUROBI_LIBS EXCLUDE REGEX "_light")

if (NOT _GUROBI_LIBS)
    message(FATAL_ERROR "GUROBI library not found in ${GUROBI_ROOT}/lib")
endif ()

list(SORT _GUROBI_LIBS)
list(LENGTH _GUROBI_LIBS _GUROBI_LIBS_LEN)
math(EXPR _GUROBI_LAST "${_GUROBI_LIBS_LEN} - 1")
list(GET _GUROBI_LIBS ${_GUROBI_LAST} GUROBI_LIBRARY)
message(STATUS "Found GUROBI library: ${GUROBI_LIBRARY}")

set(GUROBI_INCLUDE_DIRS "${GUROBI_INCLUDE_DIR}")
set(GUROBI_LIBRARIES "${GUROBI_LIBRARY}")

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GUROBI DEFAULT_MSG GUROBI_LIBRARY GUROBI_INCLUDE_DIR)


mark_as_advanced(GUROBI_INCLUDE_DIR GUROBI_LIBRARY)
