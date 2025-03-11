# FindXGB.cmake

set(XGB_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/xgb")

if (EXISTS "${XGB_ROOT}")
    set(XGB_FOUND TRUE)
    message(STATUS "Found XGB root at: ${XGB_ROOT}")
else ()
    set(XGB_FOUND FALSE)
    message(FATAL_ERROR "XGB not found at ${XGB_ROOT}")
endif ()

find_path(XGB_INCLUDE_DIR
        NAMES xgboost/c_api.h
        PATHS "${XGB_ROOT}/include"
)

find_library(XGB_LIBRARY
        NAMES libxgboost.so
        PATHS "${XGB_ROOT}/lib"
)

set(XGB_INCLUDE_DIRS "${XGB_INCLUDE_DIR}")
set(XGB_LIBRARIES "${XGB_LIBRARY}")

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(XGB DEFAULT_MSG XGB_LIBRARY XGB_INCLUDE_DIR)

mark_as_advanced(XGB_INCLUDE_DIR XGB_LIBRARY)

