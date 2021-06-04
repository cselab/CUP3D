# Look for HDF5, FFTW and GSL, all required dependencies.
# If not found, but if the CUP_AUTO_INSTALL_DEPENDENCIES is set, the
# ./install_dependencies.sh script will be run automatically.

set(_CUP_HDF5_BUILD_DIR "${DEP_BUILD_DIR}/hdf5-1.10.1")
set(_CUP_FFTW_BUILD_DIR "${DEP_BUILD_DIR}/fftw-3.3.7")
set(_CUP_GSL_BUILD_DIR "${DEP_BUILD_DIR}/gsl-2.1")

set(_missing_dep)
set(_dep_install_flags)

# Update search dirs (immediately and after running ./install_dependencies.sh).
macro(_update_dirs)
    if (NOT HDF5_DIR AND EXISTS "${_CUP_HDF5_BUILD_DIR}")
        set(HDF5_DIR "${_CUP_HDF5_BUILD_DIR}")
    endif()
    if (NOT FTTW_DIR AND NOT DEFINED ENV{FFTWDIR} AND EXISTS "${_CUP_FFTW_BUILD_DIR}")
        set(FFTW_DIR "${_CUP_FFTW_BUILD_DIR}")
    endif()
    if (NOT GSL_ROOT_DIR AND NOT DEFINED ENV{GSL_ROOT_DIR} AND EXISTS "${_CUP_GSL_BUILD_DIR}")
        set(GSL_ROOT_DIR "${_CUP_GSL_BUILD_DIR}")
    endif()
endmacro()

# Helper for adding a package to the list of missing dependency if it is not found.
macro(_find_package package flag)
    find_package(${package})
    if (NOT ${package}_FOUND)
        list(APPEND _missing_dep ${package})
        list(APPEND _dep_install_flags "${flag}")
    endif()
endmacro()

_update_dirs()

set(HDF5_PREFER_PARALLEL ON)
_find_package(HDF5 "--hdf5")

# TODO: FFTW single precision
# set(FFTW_USE_STATIC_LIBS 1)
_find_package(FFTW "--fftw")

_find_package(GSL "--gsl")

if (_missing_dep)
    string(JOIN " " _dep_install_flags_str ${_dep_install_flags})
    if (NOT CUP_AUTO_INSTALL_DEPENDENCIES)
        string(JOIN ", " _missing_dep ${_missing_dep})
        message(FATAL_ERROR
                "One or more dependencies not found: ${_missing_dep}. "
                "Dependencies can be installed by running\n"
                "    ./install_dependencies.sh ${_dep_install_flags_str}\n"
                "from the repository root, or by rerunning cmake with "
                "-DCUP_AUTO_INSTALL_DEPENDENCIES=ON")
    endif()
    message("Installing dependencies: ${_dep_install_flags_str}")
    execute_process(
        COMMAND "./install_dependencies.sh" ${_dep_install_flags}
        WORKING_DIRECTORY "${ROOT_DIR}"
        ERROR_VARIABLE error
        RESULT_VARIABLE error_code)
    if (error_code)
        message(FATAL_ERROR "Installing dependencies ${_dep_install_flags_str} failed:\n${error}")
    else()
        message("./install_dependencies.sh ${_dep_install_flags_str} completed successfully.")
    endif()
    _update_dirs()
endif()

if (NOT HDF5_FOUND)
    find_package(HDF5 REQUIRED)
endif()
if (NOT FFTW_FOUND)
    find_package(FFTW REQUIRED)
endif()
if (NOT GSL_FOUND)
    find_package(GSL REQUIRED)
endif()
