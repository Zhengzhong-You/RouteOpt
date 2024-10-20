
find_path(CVRPSEP_INCLUDE_DIR
        NAMES basegrph.h
        binpack.h
        blocks.h
        brnching.h
        capsep.h
        cnstrmgr.h
        combsep.h
        compcuts.h
        compress.h
        cutbase.h
        fcapfix.h
        fcisep.h
        fcits.h
        glmsep.h
        grsearch.h
        hpmstar.h
        htoursep.h
        intap.h
        memmod.h
        mstarsep.h
        mxf.h
        newhtour.h
        sort.h
        strcomb.h
        strngcmp.h
        tolerances.h
        twomatch.h
        PATHS "$ENV{CVRPSEP_HOME}"
        )

find_library(CVRPSEP_LIBRARY
        NAMES libcvrpsep.a
        PATHS "$ENV{CVRPSEP_HOME}/obj"
        )


set(CVRPSEP_INCLUDE_DIRS "${CVRPSEP_INCLUDE_DIR}")
set(CVRPSEP_LIBRARIES "${CVRPSEP_LIBRARY}")

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CVRPSEP DEFAULT_MSG CVRPSEP_LIBRARY CVRPSEP_INCLUDE_DIR)

mark_as_advanced(CVRPSEP_INCLUDE_DIR CVRPSEP_LIBRARY)
