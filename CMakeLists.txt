cmake_minimum_required(VERSION 3.2)

project(CubismUP_3D)

include(ProcessorCount)
ProcessorCount(NUM_CORES)

# Choose Release mode as default.
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING
        "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel."
        FORCE)
endif(NOT CMAKE_BUILD_TYPE)

# Options
set(EXENAME "Swim3D")
set(BLOCK_SIZE "16" CACHE String "Number of grid points in a block, per dimension")
set(NUM_THREADS ${NUM_CORES} CACHE String "Number of threads per node")

option(COMPILE_STATIC_LIB "Compile static library" OFF)
option(_UNBOUNDED_FFT_ "Use FFT solver for unbounded domain" OFF)
option(_DUMPGRID_ "Use asynchronous data dumps (requires additional memory)" ON)
option(_DUMP_RAW_ "Dump binary surface data for each obstacle" OFF)
option(_FLOAT_PRECISION_ "Compute in single precision" OFF)  # Because cmake/FindFFTW.cmake handles now only double precision.
option(_HDF5_DOUBLE_PRECISION_ "Dump HDF5 data in double precision" OFF)
option(_RK2_ "Use second order Runge-Kutta integrator" OFF)
option(_VERBOSE_ "Be verbose" OFF)

# Pybind11.
# Not using `python3-config` because it doesn't seem to work with `pip install --user`.
execute_process(
    COMMAND python3 -m pybind11 --includes
    OUTPUT_VARIABLE PYBIND11_INCLUDES
    OUTPUT_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE retcode1)
execute_process(
    COMMAND python3 -c "import sysconfig; print(sysconfig.get_config_var('EXT_SUFFIX')[:-3])"
    OUTPUT_VARIABLE PY_SO_SUFFIX
    OUTPUT_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE retcode2)
separate_arguments(PYBIND11_INCLUDES)
set(PYBIND11_FLAGS -shared -fPIC ${PYBIND11_INCLUDES})
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(PYBIND11_FLAGS ${PYBIND11_FLAGS} -undefined dynamic_lookup)
endif()

set(COMPILE_PY_SO "1")
if (NOT "${retcode1}" STREQUAL "0")
    message(WARNING "pybind11 not installed, not compiling .so!")
    set(COMPILE_PY_SO "0")
endif()
if (NOT "${retcode2}" STREQUAL "0")
    message(WARNING "Cannot get Python .so suffix, not compiling it!")
    set(COMPILE_PY_SO "0")
endif()


# Source files and target
set(EXE1 "${EXENAME}_open")
set(EXE2 "${EXENAME}_periodicz")
set(STATIC_LIB "${EXENAME}")
set(PY_SO "_cubismup3d${PY_SO_SUFFIX}")  # Has to match the name in `PYBIND11_MODULE`.
                                         # Wrapped by `python/cubismup3d/__init__.py`.
set(CMAKE_SHARED_LIBRARY_PREFIX "")      # No "lib" prefix for shared.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ".")
set(ROOT_FOLDER ${CMAKE_CURRENT_LIST_DIR})

include_directories(${ROOT_FOLDER}/Cubism/include/)
include_directories(${ROOT_FOLDER}/source/)
set(COMMON_SOURCES           # Common for app and lib.
    ${ROOT_FOLDER}/source/BufferedLogger.cpp
    ${ROOT_FOLDER}/source/IF3D_CarlingFishOperator.cpp
    ${ROOT_FOLDER}/source/IF3D_DCylinderObstacleOperator.cpp
    ${ROOT_FOLDER}/source/IF3D_ExternalObstacleOperator.cpp
    ${ROOT_FOLDER}/source/IF3D_FishLibrary.cpp
    ${ROOT_FOLDER}/source/IF3D_FishOperator.cpp
    ${ROOT_FOLDER}/source/IF3D_NacaOperator.cpp
    ${ROOT_FOLDER}/source/IF3D_ElasticFishOperator.cpp
    ${ROOT_FOLDER}/source/IF3D_ObstacleFactory.cpp
    ${ROOT_FOLDER}/source/IF3D_ObstacleOperator.cpp
    ${ROOT_FOLDER}/source/IF3D_ObstacleVector.cpp
    ${ROOT_FOLDER}/source/IF3D_PlateObstacleOperator.cpp
    ${ROOT_FOLDER}/source/IF3D_SphereObstacleOperator.cpp
    ${ROOT_FOLDER}/source/IF3D_StefanFishOperator.cpp
    ${ROOT_FOLDER}/source/Simulation.cpp
)

# GCC flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-parameter -Wno-vla -Wno-reorder -Wno-maybe-uninitialized")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DNDEBUG -O3 -march=native -fstrict-aliasing -ffast-math -falign-functions=32")
set(CMAKE_CXX_STANDARD 11)

# Define targets (executables). We have two: one for open, one for periodic
# boundary conditions (periodic in Z).
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${ROOT_FOLDER}/bin)
add_executable(${EXE1} ${COMMON_SOURCES} ${ROOT_FOLDER}/source/main.cpp)
add_executable(${EXE2} ${COMMON_SOURCES} ${ROOT_FOLDER}/source/main.cpp)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${ROOT_FOLDER}/lib)
if (COMPILE_STATIC_LIB)
    add_library(${STATIC_LIB}
                ${COMMON_SOURCES}
                ${ROOT_FOLDER}/source/cubism_main.cpp
                ${ROOT_FOLDER}/source/CubismUP_3D.cpp)
endif()
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${ROOT_FOLDER}/python/cubismup3d)
if ("${COMPILE_PY_SO}" STREQUAL "1")
    add_library(${PY_SO} SHARED
                ${COMMON_SOURCES}
                ${ROOT_FOLDER}/source/bindings/Simulation.cpp)
endif()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# Target-specific flags
target_compile_definitions(${EXE2} PRIVATE "-DBC_PERIODICZ")
target_compile_options(${PY_SO} PRIVATE ${PYBIND11_FLAGS})
# https://github.com/conda-forge/hdf5-feedstock/issues/58#issuecomment-310865909
target_compile_definitions(${PY_SO} PRIVATE "-DH5_BUILT_AS_DYNAMIC_LIB")

# Global flags
if (_UNBOUNDED_FFT_)
    add_definitions("-D_UNBOUNDED_FFT_")
endif()
if (_DUMP_RAW_)
    add_definitions("-D_DUMP_RAW_")
endif()
if (_FLOAT_PRECISION_)
    add_definitions("-D_FLOAT_PRECISION_")
endif()
if (_HDF5_DOUBLE_PRECISION_)
    add_definitions("-D_HDF5_DOUBLE_PRECISION_")
endif()
if (_RK2_)
    add_definitions("-D_RK2_")
endif()
if (_VERBOSE_)
    add_definitions("-D_VERBOSE_")
endif()
add_definitions("-D_ALIGNBYTES_=32")
add_definitions("-DCUBISMUP3D_BLOCK_SIZE=${BLOCK_SIZE}")
add_definitions("-DNTHREADS=${NUM_THREADS}")

# OpenMP support
find_package(OpenMP REQUIRED)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

# MPI
find_package(MPI REQUIRED)
add_definitions("-DMPI_ENABLE")
set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})

# Check if MPI thread safe
execute_process(COMMAND ${CMAKE_CXX_COMPILER} ${ROOT_FOLDER}/cmake/check_mpi_thread_safety.cpp -o check_mpi_thread_safety
                RESULT_VARIABLE _IS_COMPILED)
if (NOT _IS_COMPILED STREQUAL "0")
    message(FATAL_ERROR "Cannot compile a simple MPI program to determine whether MPI is thread safe.")
endif()
execute_process(COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ./check_mpi_thread_safety
                ERROR_QUIET
                OUTPUT_QUIET
                RESULT_VARIABLE CHECK_MPI_THREAD_SAFETY)
if (${CHECK_MPI_THREAD_SAFETY} EQUAL 0)
    if (_DUMPGRID_)
        add_definitions("-DDUMPGRID")
    endif()
else()
    message("Note: MPI not thread safe, fallback to synchronous HDF5 dumping!")
endif()


# CMake for some reason decides to replace `/full/path/libasdf.so` with
# `-lasdf` sometimes, and then messes up everything. This hack prevents that.
set(LIBPREFIX "///")

# HDF5
if (NOT HDF5_ROOT AND EXISTS ${ROOT_FOLDER}/dependencies/build/hdf5-1.10.1)
    set(HDF5_ROOT ${ROOT_FOLDER}/dependencies/build/hdf5-1.10.1)
endif()
set(HDF5_PREFER_PARALLEL ON)
find_package(HDF5 REQUIRED)
add_definitions("-DCUBISM_USE_HDF")
include_directories(${HDF5_INCLUDE_DIRS})
target_link_libraries(${EXE1} ${LIBPREFIX}${HDF5_LIBRARIES})
target_link_libraries(${EXE2} ${LIBPREFIX}${HDF5_LIBRARIES})
target_link_libraries(${PY_SO} ${LIBPREFIX}${HDF5_LIBRARIES})

# FFTW
# TODO: Single precision
if (NOT FTTW_ROOT AND NOT DEFINED ENV{FFTWDIR} AND EXISTS ${ROOT_FOLDER}/dependencies/build/fftw-3.3.7)
    set(FFTW_ROOT ${ROOT_FOLDER}/dependencies/build/fftw-3.3.7)
endif()
find_package(FFTW REQUIRED)
include_directories(${FFTW_INCLUDES})
target_link_libraries(${EXE1} ${LIBPREFIX}${FFTW_LIBRARIES})
target_link_libraries(${EXE2} ${LIBPREFIX}${FFTW_LIBRARIES})
target_link_libraries(${PY_SO} ${LIBPREFIX}${FFTW_LIBRARIES})

# GSL
if (NOT GSL_ROOT_DIR AND NOT DEFINED ENV{GSL_ROOT_DIR} AND EXISTS ${ROOT_FOLDER}/dependencies/build/gsl-2.1)
    set(GSL_ROOT_DIR ${ROOT_FOLDER}/dependencies/build/gsl-2.1)
endif()
find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})
target_link_libraries(${EXE1} ${LIBPREFIX}${GSL_LIBRARIES})
target_link_libraries(${EXE2} ${LIBPREFIX}${GSL_LIBRARIES})
target_link_libraries(${PY_SO} ${LIBPREFIX}${GSL_LIBRARIES})
