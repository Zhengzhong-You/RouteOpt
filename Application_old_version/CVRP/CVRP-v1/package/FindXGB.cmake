
find_path(XGB_INCLUDE_DIR
NAMES xgboost/c_api.h
PATHS "$ENV{XGB_HOME}/include"
)

find_library(XGB_LIBRARY
        NAMES libxgboost.so
        PATHS "$ENV{XGB_HOME}/lib"
        )


set(XGB_INCLUDE_DIRS "${XGB_INCLUDE_DIR}" )
set(XGB_LIBRARIES "${XGB_LIBRARY}" )

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(XGB DEFAULT_MSG XGB_LIBRARY XGB_INCLUDE_DIR)

#mark_as_advanced(XGB_INCLUDE_DIR XGB_LIBRARY)
