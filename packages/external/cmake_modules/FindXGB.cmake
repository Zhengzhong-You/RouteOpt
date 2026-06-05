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
        NAMES xgboost objxgboost
        PATHS
        "${XGB_ROOT}/lib"
        "${XGB_ROOT}/build/lib"
        "${XGB_ROOT}/build"
        "${XGB_ROOT}/build/Release"
        "${XGB_ROOT}/build/src/Release"
        "${XGB_ROOT}/build/src/objxgboost.dir/Release"
)

find_library(DMLC_LIBRARY
        NAMES dmlc
        PATHS
        "${XGB_ROOT}/build/dmlc-core/Release"
        "${XGB_ROOT}/build/dmlc-core"
        "${XGB_ROOT}/lib"
)

set(XGB_INCLUDE_DIRS "${XGB_INCLUDE_DIR}" "${XGB_ROOT}/rabit/include" "${XGB_ROOT}/dmlc-core/include")
if (DMLC_LIBRARY)
    set(XGB_LIBRARIES "${XGB_LIBRARY}" "${DMLC_LIBRARY}")
else()
    set(XGB_LIBRARIES "${XGB_LIBRARY}")
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(XGB DEFAULT_MSG XGB_LIBRARY XGB_INCLUDE_DIR)

mark_as_advanced(XGB_INCLUDE_DIR XGB_LIBRARY DMLC_LIBRARY)

