find_path(GUROBI_INCLUDE_DIR
    NAMES gurobi_c++.h gurobi_c.h
    PATHS
        "$ENV{GUROBI_HOME}/include")

find_library( GUROBI_LIBRARY
        NAMES libgurobi100.so
        PATHS "$ENV{GUROBI_HOME}/lib"
        )

find_library( GUROBI_CXX_LIBRARY
        NAMES libgurobi_c++.a
        PATHS "$ENV{GUROBI_HOME}/lib"
        )

set(GUROBI_INCLUDE_DIRS "${GUROBI_INCLUDE_DIR}" )
set(GUROBI_LIBRARIES "${GUROBI_LIBRARY};${GUROBI_CXX_LIBRARY}" )

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GUROBI DEFAULT_MSG GUROBI_LIBRARY GUROBI_CXX_LIBRARY GUROBI_INCLUDE_DIR)

mark_as_advanced(GUROBI_INCLUDE_DIR GUROBI_LIBRARY GUROBI_CXX_LIBRARY)