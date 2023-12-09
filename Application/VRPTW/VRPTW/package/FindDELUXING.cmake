
find_path(DELUXING_INCLUDE_DIR
        NAMES deLuxing.hpp
        PATHS "$ENV{DELUXING_HOME}/include"
        )

find_library(DELUXING_LIBRARY
        NAMES libdeluxing.a
        PATHS "$ENV{DELUXING_HOME}/lib"
        )


set(DELUXING_INCLUDE_DIRS "${DELUXING_INCLUDE_DIR}")
set(DELUXING_LIBRARIES "${DELUXING_LIBRARY}")

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(DELUXING DEFAULT_MSG DELUXING_LIBRARY DELUXING_INCLUDE_DIR)

#mark_as_advanced(DELUXING_INCLUDE_DIR DELUXING_LIBRARY)
