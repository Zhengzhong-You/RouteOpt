find_path(HGS_INCLUDE_DIR
        NAMES
        my_hgs.hpp
        PATHS "$ENV{HGS_HOME}/include"
        )

find_library(HGS_LIBRARY
        NAMES libhgscvrp-yzz.so
        PATHS "$ENV{HGS_HOME}/lib"
        )


set(HGS_INCLUDE_DIRS "${HGS_INCLUDE_DIR}")
set(HGS_LIBRARIES "${HGS_LIBRARY}")

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(HGS DEFAULT_MSG HGS_LIBRARY HGS_INCLUDE_DIR)

#mark_as_advanced(HGS_INCLUDE_DIR HGS_LIBRARY)
