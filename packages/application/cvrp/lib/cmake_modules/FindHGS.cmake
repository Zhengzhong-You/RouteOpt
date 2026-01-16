set(HGS_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/lib/hgs")

find_path(HGS_INCLUDE_DIR
        NAMES
        my_hgs.hpp
        PATHS "${HGS_ROOT}/include"
)

find_library(HGS_LIBRARY
        NAMES hgscvrp
        PATHS
        "${HGS_ROOT}/lib"
        "${HGS_ROOT}/build/lib"
)


set(HGS_INCLUDE_DIRS "${HGS_INCLUDE_DIR}")
set(HGS_LIBRARIES "${HGS_LIBRARY}")

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(HGS DEFAULT_MSG HGS_LIBRARY HGS_INCLUDE_DIR)

#mark_as_advanced(HGS_INCLUDE_DIR HGS_LIBRARY)
