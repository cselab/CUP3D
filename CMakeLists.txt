cmake_minimum_required(VERSION 3.2)

project(CubismUP_3D)

include(ProcessorCount)
ProcessorCount(NUM_CORES)

# Choose Release mode as default.
if (NOT EXISTS ${CMAKE_BINARY_DIR}/CMakeCache.txt)
  if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
  endif()
endif()

# Options
set(BLOCK_SIZE "16" CACHE String "Number of grid points in a block, per dimension")
set(NUM_THREADS ${NUM_CORES} CACHE String "Number of threads per node")
# option(SINGLE_PRECISION OFF)
set(SINGLE_PRECISION OFF)  # Because cmake/FindFFTW.cmake handles now only double precision.
option(VERBOSE OFF)

# Source files and target
set(EXE1 "cubismup3d_open")
set(EXE2 "cubismup3d_periodicz")
set(LIB "cubismup3d")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ".")
set(ROOT_FOLDER ${CMAKE_CURRENT_LIST_DIR})

include_directories(${ROOT_FOLDER}/Cubism/source/)
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
    ${ROOT_FOLDER}/source/IF3D_VortexOperator.cpp
    ${ROOT_FOLDER}/source/Simulation.cpp
)

# Define targets (executables). We have two: one for open, one for periodic
# boundary conditions (periodic in Z).
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${ROOT_FOLDER}/bin)
add_executable(${EXE1} ${COMMON_SOURCES} ${ROOT_FOLDER}/source/main.cpp)
add_executable(${EXE2} ${COMMON_SOURCES} ${ROOT_FOLDER}/source/main.cpp)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${ROOT_FOLDER}/lib)
add_library(${LIB}
            ${COMMON_SOURCES}
            ${ROOT_FOLDER}/source/cubism_main.cpp
            ${ROOT_FOLDER}/source/CubismUP_3D.cpp)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# Executable-specific flags
target_compile_definitions(${EXE2} PRIVATE "-DBC_PERIODICZ")

# Global flags
if (SINGLE_PRECISION)
    add_definitions("-D_FLOAT_PRECISION_")
endif()
if (VERBOSE)
    add_definitions("-D_VERBOSE_")
endif()
add_definitions("-D_ALIGNBYTES_=32")
add_definitions("-D_BLOCKSIZE_=${BLOCK_SIZE}")
add_definitions("-D_BS_=${BLOCK_SIZE}")
add_definitions("-D_BLOCKSIZEX_=${BLOCK_SIZE}")
add_definitions("-D_BLOCKSIZEY_=${BLOCK_SIZE}")
add_definitions("-D_BLOCKSIZEZ_=${BLOCK_SIZE}")
add_definitions("-DNTHREADS=${NUM_THREADS}")

# GCC flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-parameter -Wno-vla -Wno-reorder -Wno-maybe-uninitialized")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -march=native -DNDEBUG -fstrict-aliasing -ffast-math -falign-functions=32")
set(CMAKE_CXX_STANDARD 11)

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
    add_definitions("-DDUMPGRID")
else()
    message("Note: MPI not thread safe, fallback to synchronous HDF5 dumping!")
endif()


# HDF5
set(HDF5_PREFER_PARALLEL ON)
find_package(HDF5 REQUIRED)
add_definitions("-D_USE_HDF_")
include_directories(${HDF5_INCLUDE_DIRS})
target_link_libraries(${EXE1} ${HDF5_LIBRARIES})
target_link_libraries(${EXE2} ${HDF5_LIBRARIES})

# FFTW
# TODO: Single precision
find_package(FFTW REQUIRED)
include_directories(${FFTW_INCLUDES})
target_link_libraries(${EXE1} ${FFTW_LIBRARIES})
target_link_libraries(${EXE2} ${FFTW_LIBRARIES})

# GSL
find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})
target_link_libraries(${EXE1} ${GSL_LIBRARIES})
target_link_libraries(${EXE2} ${GSL_LIBRARIES})
