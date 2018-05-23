cmake_minimum_required(VERSION 3.2)

project(CubismUP_3D)

# Options
# option(SINGLE_PRECISION OFF)
set(SINGLE_PRECISION OFF)  # Because cmake/FindFFTW.cmake handles now only double precision.

option(VERBOSE OFF)
set(BOUNDARY_CONDITIONS "open" CACHE String "Either 'open' or 'periodicz'")
set(BLOCK_SIZE "16" CACHE String "Number of grid points in a block, per dimension")
set(NUM_THREADS "24" CACHE String "Number of threads per node")

# Source files and target
set(EXE "cubismup_3d")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "makefiles/")

include_directories(Cubism/)
set(COMMON_SOURCES           # Common for app and lib.
    # source/IF3D_CarlingFishOperator.cpp
    source/IF3D_DCylinderObstacleOperator.cpp
    # source/IF3D_DeadFishOperator.cpp
    source/IF3D_FishLibrary.cpp
    source/IF3D_FishOperator.cpp
    source/IF3D_ForcedSphereObstacleOperator.cpp
    source/IF3D_NacaOperator.cpp
    source/IF3D_ObstacleFactory.cpp
    source/IF3D_ObstacleOperator.cpp
    source/IF3D_ObstacleVector.cpp
    source/IF3D_PlateObstacleOperator.cpp
    source/IF3D_SphereObstacleOperator.cpp
    source/IF3D_StefanFishOperator.cpp
    source/IF3D_VortexOperator.cpp
    source/ProcessOperatorsOMP.cpp
    source/Save_splicer.cpp
    source/Simulation.cpp
    Cubism/Profiler.cpp

    # Only if not ACCFFT.
    source/PoissonSolverScalarFFTW.cpp
    source/PoissonSolverScalarFFTW_MPI.cpp
)
add_executable(${EXE}
               ${COMMON_SOURCES}
               source/main.cpp)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# Project-specific flags
if (SINGLE_PRECISION)
    add_definitions("-D_SP_COMP_")
endif()
if (VERBOSE)
    add_definitions("-D_VERBOSE_")
endif()
add_definitions("-D_ALIGNBYTES_=32")
add_definitions("-D_BLOCKSIZE_=${BLOCK_SIZE}")
add_definitions("-D_BS_=${BLOCK_SIZE}")
add_definitions("-D_BSX_=${BLOCK_SIZE}")
add_definitions("-D_BSY_=${BLOCK_SIZE}")
add_definitions("-D_BSZ_=${BLOCK_SIZE}")
add_definitions("-DNTHREADS=${NUM_THREADS}")

# GCC flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-parameter -Wno-vla -Wno-reorder -Wno-maybe-uninitialized")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -march=native -DNDEBUG -fstrict-aliasing -ffast-math -falign-functions=32")

set(CMAKE_CXX_STANDARD 11)

# OpenMP support
find_package(OpenMP)
if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS
        "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# MPI
find_package(MPI REQUIRED)
add_definitions("-DMPI_ENABLE")
set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})

# HDF5
set(HDF5_PREFER_PARALLEL ON)
find_package(HDF5 REQUIRED)
add_definitions("-D_USE_HDF_")
include_directories(${HDF5_INCLUDE_DIRS})
target_link_libraries(${EXE} ${HDF5_LIBRARIES})

# FFTW
# TODO: Single precision
find_package(FFTW REQUIRED)
message("Found: ${FFTW_FOUND} include: ${FFTW_INCLUDES} libs: ${FFTW_LIBRARIES}")
include_directories(${FFTW_INCLUDES})
target_link_libraries(${EXE} ${FFTW_LIBRARIES})

# GSL
find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDES})
target_link_libraries(${EXE} ${GSL_LIBRARIES})
