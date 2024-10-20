

find_path(CPLEX_INCLUDE_DIR
        NAMES ilcplex/cplex.h
        PATHS "$ENV{CPLEX_HOME}/cplex/include"
        #        PATH_SUFFIXES include
        NO_DEFAULT_PATH
)

find_library(CPLEX_LIBRARY
        NAMES cplex
        PATHS "$ENV{CPLEX_HOME}/cplex/lib"
        PATH_SUFFIXES x86-64_linux/static_pic x64_windows/static_md x86_windows/static_md arm64_osx/static_pic
        NO_DEFAULT_PATH
)

find_library(CPLEX_CONCERT_LIBRARY
        NAMES concert
        PATHS "$ENV{CPLEX_HOME}/concert/lib"
        PATH_SUFFIXES x86-64_linux/static_pic x64_windows/static_md x86_windows/static_md arm64_osx/static_pic
        NO_DEFAULT_PATH
)

find_library(CPLEX_ILOCPLEX_LIBRARY
        NAMES ilocplex
        PATHS "$ENV{CPLEX_HOME}/cplex/lib"
        PATH_SUFFIXES x86-64_linux/static_pic x64_windows/static_md x86_windows/static_md arm64_osx/static_pic
        NO_DEFAULT_PATH
)

set(CPLEX_INCLUDE_DIRS ${CPLEX_INCLUDE_DIR})
set(CPLEX_LIBRARIES ${CPLEX_LIBRARY} ${CPLEX_CONCERT_LIBRARY} ${CPLEX_ILOCPLEX_LIBRARY})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CPLEX DEFAULT_MSG CPLEX_LIBRARY CPLEX_CONCERT_LIBRARY CPLEX_ILOCPLEX_LIBRARY CPLEX_INCLUDE_DIR)

mark_as_advanced(CPLEX_INCLUDE_DIR CPLEX_LIBRARY CPLEX_CONCERT_LIBRARY CPLEX_ILOCPLEX_LIBRARY)
