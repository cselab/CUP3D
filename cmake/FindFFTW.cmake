# https://github.com/libigl/eigen/blob/master/cmake/FindFFTW.cmake
# (modified)

# - Find the FFTW library
#
# Usage:
#   find_package(FFTW [REQUIRED] [QUIET] )
#
# It sets the following variables:
#   FFTW_FOUND               ... true if fftw is found on the system
#   FFTW_LIBRARIES           ... full path to fftw library
#   FFTW_INCLUDES            ... fftw include directory
#
# The following variables will be checked by the function
#   FFTW_USE_STATIC_LIBS    ... if true, only static libraries are found
#   FFTW_ROOT               ... if set, the libraries are exclusively searched
#                               under this path
#   FFTW_LIBRARY            ... fftw library to use
#   FFTW_INCLUDE_DIR        ... fftw include directory
#

# If environment variable FFTWDIR is specified, it has the same effect as FFTW_ROOT.
if (NOT FFTW_ROOT AND DEFINED ENV{FFTWDIR})
    set( FFTW_ROOT $ENV{FFTWDIR} )
endif()

# Check if we can use PkgConfig
find_package(PkgConfig)

#Check whether to search static or dynamic libs
set(CMAKE_FIND_LIBRARY_SUFFIXES_SAV ${CMAKE_FIND_LIBRARY_SUFFIXES})

if (${FFTW_USE_STATIC_LIBS})
    set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_STATIC_LIBRARY_SUFFIX})
else()
    set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_SHARED_LIBRARY_SUFFIX})
endif()

if (FFTW_ROOT)
    find_library(FFTW_LIB1 NAMES "fftw3" PATHS ${FFTW_ROOT} PATH_SUFFIXES "lib" NO_DEFAULT_PATH)
    find_library(FFTW_LIB2 NAMES "fftw3_threads" PATHS ${FFTW_ROOT} PATH_SUFFIXES "lib" NO_DEFAULT_PATH)
    find_library(FFTW_LIB3 NAMES "fftw3_mpi" PATHS ${FFTW_ROOT} PATH_SUFFIXES "lib" NO_DEFAULT_PATH)
    find_path(FFTW_INCLUDES NAMES "fftw3-mpi.h" PATHS ${FFTW_ROOT} PATH_SUFFIXES "include" NO_DEFAULT_PATH)
else()
# Determine from PKG.
    if (PKG_CONFIG_FOUND AND NOT FFTW_ROOT)
        pkg_check_modules(PKG_FFTW QUIET "fftw3")
    endif()

    find_library(FFTW_LIB1 NAMES "fftw3" PATHS ${PKG_FFTW_LIBRARY_DIRS} ${LIB_INSTALL_DIR})
    find_library(FFTW_LIB2 NAMES "fftw3_threads" PATHS ${PKG_FFTW_LIBRARY_DIRS} ${LIB_INSTALL_DIR})
    find_library(FFTW_LIBr3 NAMES "fftw3_mpi" PATHS ${PKG_FFTW_LIBRARY_DIRS} ${LIB_INSTALL_DIR})
    find_path(
        FFTW_INCLUDES
        NAMES "fftw3-mpi.h"
        PATHS ${PKG_FFTW_INCLUDE_DIRS} ${INCLUDE_INSTALL_DIR}
    )
endif(FFTW_ROOT)

set(FFTW_LIBRARIES ${FFTW_LIB1} ${FFTW_LIB2} ${FFTW_LIB3})

set( CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES_SAV} )

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FFTW DEFAULT_MSG
                                  FFTW_INCLUDES FFTW_LIBRARIES)

mark_as_advanced(FFTW_INCLUDES FFTW_LIBRARIES FFTW_LIB FFTWF_LIB FFTWL_LIB)

