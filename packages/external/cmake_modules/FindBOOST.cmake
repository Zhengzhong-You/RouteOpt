# FindBOOST.cmake

set(BOOST_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/boost")

if (EXISTS "${BOOST_ROOT}")
    set(BOOST_FOUND TRUE)
    message(STATUS "Found BOOST root at: ${BOOST_ROOT}")
else ()
    set(BOOST_FOUND FALSE)
    message(FATAL_ERROR "BOOST not found at ${BOOST_ROOT}")
endif ()

# Directly include all headers under BOOST/
set(BOOST_INCLUDE_DIRS "${BOOST_ROOT}")


include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(BOOST DEFAULT_MSG BOOST_INCLUDE_DIRS)

mark_as_advanced(BOOST_INCLUDE_DIRS)
