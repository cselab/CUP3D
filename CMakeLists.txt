# NOTE: Update `install_dependencies.sh` if updating the minimum required version.
# Minimum version 3.9 for find_package(MPI).
# Minimum version 3.12 for find_package(Python).
cmake_minimum_required(VERSION 3.12)
project(CubismUP_3D)

if(POLICY CMP0074)
    cmake_policy(SET CMP0074 NEW)
endif()

# Set Release mode as default.
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING
        "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel."
        FORCE)
endif(NOT CMAKE_BUILD_TYPE)

# Paths and target names.
set(BUILD_DIR "${CMAKE_CURRENT_BINARY_DIR}")
set(ROOT_DIR "${CMAKE_CURRENT_LIST_DIR}")
set(GENERATED_DIR "${BUILD_DIR}/_generated")
set(SRC_DIR "${ROOT_DIR}/source")
set(DEP_BUILD_DIR "${ROOT_DIR}/dependencies/build")

set(EXE "cubismup3d_simulation")
set(CORE "cubismup3d_core")
set(CORE_PIC "cubismup3d_core_pic")  # Position independent code.
set(PUBLICDEP "cubismup3d_publicdep")    # Public dependencies and flags.
set(PRIVATEDEP "cubismup3d_privatedep")  # Private dependencies and flags.
set(LIB "cubismup3d")
set(PYLIB "libcubismup3d")


# Options.
set(CUP_BLOCK_SIZE "16" CACHE STRING "Number of grid points in a block, per dimension")
option(CUP_TESTS "Enable tests" OFF)
option(CUP_BUILD_EXE "Build the executable" ON)
option(CUP_BUILD_LIB "Build the static library" OFF)
option(CUP_BUILD_PY "Build Python bindings" OFF)

# Compile-time settings. Stored in a configuration file CubismUP3DMacros.h and imported by Definitions.h.
option(CUP_ASYNC_DUMP "Use asynchronous data dumps (requires additional memory)" ON)
option(CUP_DUMP_SURFACE_BINARY "Dump binary surface data for each obstacle" OFF)
option(CUP_HDF5_DOUBLE_PRECISION "Dump HDF5 data in double precision" OFF)
option(CUP_SINGLE_PRECISION "Compute in single precision" OFF)  # Because cmake/FindFFTW.cmake handles now only double precision.
option(CUP_VERBOSE "Be verbose" OFF)

if (CUP_SINGLE_PRECISION)
    message(FATAL_ERROR "Single precision compilation not implemented. Check cmake/FindFFTW.cmake first.")
endif()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# MPI - must be before modifying CMAKE_SHARED_LIBRARY_SUFFIX (?).
# Check if MPI is thread safe.
find_package(MPI REQUIRED)
execute_process(COMMAND ${MPI_CXX_COMPILER} ${ROOT_DIR}/cmake/check_mpi_thread_safety.cpp -o check_mpi_thread_safety
                RESULT_VARIABLE _IS_COMPILED)
if (NOT _IS_COMPILED STREQUAL "0")
    message(FATAL_ERROR "Cannot compile a simple MPI program to determine whether MPI is thread safe.")
endif()
execute_process(COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ./check_mpi_thread_safety
                ERROR_QUIET
                OUTPUT_QUIET
                RESULT_VARIABLE CUP_IS_MPI_THREAD_UNSAFE)
if (${CUP_IS_MPI_THREAD_UNSAFE})
    message("Note: MPI not thread safe, fallback to synchronous HDF5 dumping!")
endif()

# Cubism.
add_subdirectory("${ROOT_DIR}/Cubism")

# OpenMP.
find_package(OpenMP REQUIRED)

# HDF5.
if (NOT HDF5_DIR AND EXISTS ${ROOT_FOLDER}/dependencies/build/hdf5-1.10.1)
    set(HDF5_DIR ${ROOT_FOLDER}/dependencies/build/hdf5-1.10.1)
endif()
set(HDF5_PREFER_PARALLEL ON)
find_package(HDF5 REQUIRED)

# FFTW.
# TODO: Single precision
if (NOT FTTW_DIR AND NOT DEFINED ENV{FFTWDIR} AND EXISTS ${DEP_BUILD_DIR}/fftw-3.3.7)
    set(FFTW_DIR ${DEP_BUILD_DIR}/fftw-3.3.7)
endif()
# set(FFTW_USE_STATIC_LIBS 1)
find_package(FFTW REQUIRED)

# GSL.
if (NOT GSL_ROOT_DIR AND NOT DEFINED ENV{GSL_ROOT_DIR} AND EXISTS ${DEP_BUILD_DIR}/gsl-2.1)
    set(GSL_ROOT_DIR ${DEP_BUILD_DIR}/gsl-2.1)
endif()
find_package(GSL REQUIRED)


# Targets.
set(COMMON_SOURCES           # Common for app and lib.
    "${SRC_DIR}/obstacles/CarlingFish.cpp"
    "${SRC_DIR}/obstacles/Cylinder.cpp"
    # "${SRC_DIR}/obstacles/ElasticFish.cpp"
    "${SRC_DIR}/obstacles/Ellipsoid.cpp"
    "${SRC_DIR}/obstacles/ExternalObstacle.cpp"
    "${SRC_DIR}/obstacles/Fish.cpp"
    "${SRC_DIR}/obstacles/FishLibrary.cpp"
    "${SRC_DIR}/obstacles/FishShapes.cpp"
    "${SRC_DIR}/obstacles/Naca.cpp"
    "${SRC_DIR}/obstacles/Obstacle.cpp"
    "${SRC_DIR}/obstacles/ObstacleFactory.cpp"
    "${SRC_DIR}/obstacles/ObstacleVector.cpp"
    "${SRC_DIR}/obstacles/Plate.cpp"
    "${SRC_DIR}/obstacles/Sphere.cpp"
    "${SRC_DIR}/obstacles/StefanFish.cpp"
    "${SRC_DIR}/operators/AdvectionDiffusion.cpp"
    "${SRC_DIR}/operators/Analysis.cpp"
    "${SRC_DIR}/operators/Checkpoint.cpp"
    "${SRC_DIR}/operators/ComputeDissipation.cpp"
    "${SRC_DIR}/operators/ExternalForcing.cpp"
    "${SRC_DIR}/operators/FadeOut.cpp"
    "${SRC_DIR}/operators/FixedMassFlux_nonUniform.cpp"
    "${SRC_DIR}/operators/FluidSolidForces.cpp"
    "${SRC_DIR}/operators/HITfiltering.cpp"
    "${SRC_DIR}/operators/InitialConditions.cpp"
    "${SRC_DIR}/operators/IterativePressureNonUniform.cpp"
    "${SRC_DIR}/operators/IterativePressurePenalization.cpp"
    "${SRC_DIR}/operators/ObstaclesCreate.cpp"
    "${SRC_DIR}/operators/ObstaclesUpdate.cpp"
    "${SRC_DIR}/operators/Penalization.cpp"
    "${SRC_DIR}/operators/PressureProjection.cpp"
    "${SRC_DIR}/operators/PressureRHS.cpp"
    "${SRC_DIR}/operators/SGS.cpp"
    "${SRC_DIR}/poisson/PoissonSolver.cpp"
    "${SRC_DIR}/poisson/PoissonSolverMixed.cpp"
    "${SRC_DIR}/poisson/PoissonSolverPeriodic.cpp"
    "${SRC_DIR}/poisson/PoissonSolverUnbounded.cpp"
    "${SRC_DIR}/spectralOperators/SpectralAnalysis.cpp"
    "${SRC_DIR}/spectralOperators/SpectralForcing.cpp"
    "${SRC_DIR}/spectralOperators/SpectralIcGenerator.cpp"
    "${SRC_DIR}/spectralOperators/SpectralManip.cpp"
    "${SRC_DIR}/spectralOperators/SpectralManipFFTW.cpp"
    "${SRC_DIR}/utils/BufferedLogger.cpp"
    "${SRC_DIR}/Simulation.cpp"
    "${SRC_DIR}/SimulationData.cpp"
)

add_library(${PUBLICDEP} INTERFACE)
target_link_libraries(${PUBLICDEP} INTERFACE OpenMP::OpenMP_CXX)
target_link_libraries(${PUBLICDEP} INTERFACE MPI::MPI_CXX)
target_link_libraries(${PUBLICDEP} INTERFACE Cubism::Cubism)
target_link_libraries(${PUBLICDEP} INTERFACE FFTW::FFTW)
target_link_libraries(${PUBLICDEP} INTERFACE GSL::gsl)
target_link_libraries(${PUBLICDEP} INTERFACE ${HDF5_LIBRARIES})
target_include_directories(${PUBLICDEP} INTERFACE ${HDF5_INCLUDE_DIRS})
target_include_directories(${PUBLICDEP} INTERFACE ${GENERATED_DIR})

add_library(${PRIVATEDEP} INTERFACE)
set(_PRIVATE_FLAGS -Wall -Wextra -Wno-unused-parameter -Wno-maybe-uninitialized -fmax-errors=3)
set(_PRIVATE_RELEASE_FLAGS -march=native -fstrict-aliasing -ffast-math -falign-functions=32)
target_compile_options(${PRIVATEDEP} INTERFACE ${_PRIVATE_FLAGS} $<$<CONFIG:Release>:${_PRIVATE_RELEASE_FLAGS}>)

add_library(${CORE} STATIC EXCLUDE_FROM_ALL ${COMMON_SOURCES})
add_library(cubismup3d::core ALIAS ${CORE})
target_link_libraries(${CORE} PRIVATE ${PRIVATEDEP})
target_link_libraries(${CORE} PUBLIC ${PUBLICDEP})

# Add a copy of the core with -fPIC enabled. PIC affects performance but is
# needed for compiling shared objects. Exclude from default make.
add_library(${CORE_PIC} STATIC EXCLUDE_FROM_ALL ${COMMON_SOURCES})
add_library(cubismup3d::core_pic ALIAS ${CORE_PIC})
target_link_libraries(${CORE_PIC} PRIVATE ${PRIVATEDEP})
target_link_libraries(${CORE_PIC} PUBLIC ${PUBLICDEP})
set_property(TARGET ${CORE_PIC} PROPERTY POSITION_INDEPENDENT_CODE 1)

# set(CMAKE_CXX_EXTENSIONS OFF)  # -std=c++17, not -std=gnu++17!
# if (NOT CMAKE_CXX_STANDARD)
#     set(CMAKE_CXX_STANDARD 17)
# endif()

if (CUP_BUILD_EXE)
    add_executable(${EXE} "${ROOT_DIR}/source/main.cpp")
    target_link_libraries(${EXE} PRIVATE ${CORE})
    target_link_libraries(${EXE} PRIVATE ${PRIVATEDEP})
endif()

if (CUP_BUILD_LIB)
    add_library(${LIB} STATIC "${ROOT_DIR}/source/cubism_main.cpp")
    target_link_libraries(${LIB} PUBLIC ${CORE})
    target_link_libraries(${LIB} PRIVATE ${PRIVATEDEP})
endif()

if (CUP_BUILD_PY)
    set(PYTHON_SOURCES
        ${SRC_DIR}/bindings/Obstacles.cpp
        ${SRC_DIR}/bindings/Simulation.cpp
    )

    # Check if pybind11 was already imported by some of the parent cmake files.
    if (TARGET pybind11::pybind11)
        # If so, directly use pybind11_add_module.
        pybind11_add_module(
            ${PYLIB}
            ${PYTHON_SOURCES})
    else()
        # Otherwise, find Python & pybind11 package, and use Python_add_library.
        find_package(Python COMPONENTS Interpreter Development)
        if (NOT ${Python_Interpreter_FOUND})
            message(FATAL_ERROR "Python not found.")
        endif()
        Python_add_library(
            ${PYLIB}
            MODULE
            WITH_SOABI
            ${PYTHON_SOURCES})
        execute_process(
            COMMAND ${Python_EXECUTABLE} -m pybind11 --includes
            OUTPUT_VARIABLE _PYBIND11_INCLUDES
            OUTPUT_STRIP_TRAILING_WHITESPACE
            RESULT_VARIABLE _RETCODE)
        if (NOT "${_RETCODE}" STREQUAL "0")
            message(FATAL_ERROR "pybind11 not installed, cannot compile Python bindings!")
        endif()
        target_include_directories(${PYLIB} PRIVATE ${_PYBIND11_INCLUDES})
    endif()
    target_link_libraries(${PYLIB} PUBLIC ${CORE_PIC})
    target_link_libraries(${PYLIB} PRIVATE ${PRIVATEDEP})

    # Old flag, not sure if needed anymore.
    # # https://github.com/conda-forge/hdf5-feedstock/issues/58#issuecomment-310865909
    # target_compile_definitions(${PYLIB} PRIVATE "-DH5_BUILT_AS_DYNAMIC_LIB")
endif()

if (CUP_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

# Generate macro file with current compilation settings. This file is generated
# at configuration time (`cmake ..`), so it should be available for compilation
# immediately.
# https://cmake.org/cmake/help/latest/command/configure_file.html
configure_file(${ROOT_DIR}/cmake/CubismUP3DMacros.template.h
               ${GENERATED_DIR}/generated.CubismUP3DMacros.h)
