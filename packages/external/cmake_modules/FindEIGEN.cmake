# FindEIGEN.cmake

set(EIGEN_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/eigen")

if (EXISTS "${EIGEN_ROOT}")
    set(EIGEN_FOUND TRUE)
    message(STATUS "Found EIGEN root at: ${EIGEN_ROOT}")
else ()
    set(EIGEN_FOUND FALSE)
    message(FATAL_ERROR "EIGEN not found at ${EIGEN_ROOT}")
endif ()

# Directly include all headers under EIGEN/
set(EIGEN_INCLUDE_DIRS "${EIGEN_ROOT}")


include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(EIGEN DEFAULT_MSG EIGEN_INCLUDE_DIRS)

mark_as_advanced(EIGEN_INCLUDE_DIRS)

