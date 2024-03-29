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
set(PUBLICDEP "cubismup3d_publicdep")    # Public dependencies and flags.
set(PRIVATEDEP "cubismup3d_privatedep")  # Private dependencies and flags.
set(LIB "cubismup3d")
set(PYLIB "libcubismup3d")


# Options.
set(CUP_BLOCK_SIZE "8" CACHE STRING "Number of grid points in a block, per dimension")
set(CUP_ROOT_DIR "${ROOT_DIR}" CACHE INTERNAL "CubismUP_3D repository root path")
set(CUP_BUILD_DIR "${BUILD_DIR}" CACHE INTERNAL "CubismUP_3D build path")
option(CUP_BUILD_EXE "Build the executable" OFF)
option(CUP_BUILD_LIB "Build the static library" OFF)
option(CUP_BUILD_PY "Build Python bindings" ON)
option(CUP_TESTS "Enable tests" OFF)

# Compile-time settings. Stored in a configuration file CubismUP3DMacros.h and imported by Definitions.h.
option(CUP_ASYNC_DUMP "Use asynchronous data dumps (uses additional memory, requires multithreaded MPI)" ON)
option(CUP_DUMP_SURFACE_BINARY "Dump binary surface data for each obstacle" OFF)
option(CUP_HDF5_DOUBLE_PRECISION "Dump HDF5 data in double precision" OFF)
option(CUP_OPENMP "Build with OpenMP" ON)
option(CUP_VERBOSE "Be verbose" OFF)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# MPI - must be before modifying CMAKE_SHARED_LIBRARY_SUFFIX (?).
find_package(MPI REQUIRED)

# Cubism.
add_subdirectory("${ROOT_DIR}/Cubism")
# TODO: Add prefix to the DIMENSION macro.
target_compile_definitions(Cubism INTERFACE -DCUBISM_USE_MAP)
target_compile_definitions(Cubism PUBLIC CUBISM_USE_MAP)
target_compile_definitions(Cubism PUBLIC DIMENSION=3)

# OpenMP.
find_package(OpenMP REQUIRED)

# 3rd party dependencies.
include("${ROOT_DIR}/cmake/Dependencies.cmake")

# Targets.
set(CORE_SOURCE
    "${SRC_DIR}/obstacles/CarlingFish.cpp"
    "${SRC_DIR}/obstacles/Cylinder.cpp"
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
    "${SRC_DIR}/operators/ComputeDissipation.cpp"
    "${SRC_DIR}/operators/FluidSolidForces.cpp"
    "${SRC_DIR}/operators/InitialConditions.cpp"
    "${SRC_DIR}/operators/ObstaclesCreate.cpp"
    "${SRC_DIR}/operators/ObstaclesUpdate.cpp"
    "${SRC_DIR}/operators/Penalization.cpp"
    "${SRC_DIR}/operators/PressureProjection.cpp"
    "${SRC_DIR}/poisson/PoissonSolverAMR.cpp"
    "${SRC_DIR}/poisson/PoissonSolverAMRKernels.cpp"
    "${SRC_DIR}/utils/BufferedLogger.cpp"
    "${SRC_DIR}/Simulation.cpp"
    "${SRC_DIR}/SimulationData.cpp"
)

add_library(${PUBLICDEP} INTERFACE)
target_compile_definitions(${PUBLICDEP} INTERFACE -D_DOUBLE_PRECISION_)
target_link_libraries(${PUBLICDEP} INTERFACE MPI::MPI_CXX)
target_link_libraries(${PUBLICDEP} INTERFACE Cubism::Cubism)
target_link_libraries(${PUBLICDEP} INTERFACE GSL::gsl)
target_link_libraries(${PUBLICDEP} INTERFACE ${HDF5_LIBRARIES})
target_include_directories(${PUBLICDEP} INTERFACE ${HDF5_INCLUDE_DIRS})
target_include_directories(${PUBLICDEP} INTERFACE ${GENERATED_DIR})
if (CUP_OPENMP)
    target_link_libraries(${PUBLICDEP} INTERFACE OpenMP::OpenMP_CXX)
    set(_PRIVATE_FLAGS "")
else()
    target_link_libraries(${PUBLICDEP} INTERFACE ${OpenMP_CXX_LIBRARIES})
    set(_PRIVATE_FLAGS -Wno-unknown-pragmas)
endif()

add_library(${PRIVATEDEP} INTERFACE)
set(_PRIVATE_FLAGS -Wall -Wextra -Wno-unused-parameter -Wno-maybe-uninitialized -fmax-errors=3 ${_PRIVATE_FLAGS})
set(_PRIVATE_RELEASE_FLAGS
    -march=native
    -fstrict-aliasing
    -ffast-math
    -falign-functions=32)
target_compile_features(${PRIVATEDEP} INTERFACE cxx_std_17)
target_compile_options(${PRIVATEDEP} INTERFACE
                       ${_PRIVATE_FLAGS}
                       $<$<CONFIG:RelWithDebInfo>:${_PRIVATE_RELEASE_FLAGS}>
                       $<$<CONFIG:Release>:${_PRIVATE_RELEASE_FLAGS}>)

add_library(${CORE} STATIC EXCLUDE_FROM_ALL ${CORE_SOURCE})
add_library(cubismup3d::core ALIAS ${CORE})
set_property(TARGET ${CORE} PROPERTY CXX_EXTENSIONS OFF)
# PIC affects performance but is needed for compiling shared objects.
set_property(TARGET ${CORE} PROPERTY POSITION_INDEPENDENT_CODE 1)
target_link_libraries(${CORE} PRIVATE ${PRIVATEDEP})
target_link_libraries(${CORE} PUBLIC ${PUBLICDEP})
# target_compile_options(${CORE} PUBLIC -fsanitize=address -static-libasan)
# target_link_options(${CORE} PUBLIC -fsanitize=address -static-libasan)

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
    if (NOT TARGET pybind11::pybind11)
        find_package(Python COMPONENTS REQUIRED Interpreter Development)
        execute_process(
            COMMAND ${Python_EXECUTABLE} -m pybind11 --cmakedir
            OUTPUT_VARIABLE _PYBIND11_CMAKEDIR
            OUTPUT_STRIP_TRAILING_WHITESPACE
            RESULT_VARIABLE _RETCODE)
        if (NOT "${_RETCODE}" STREQUAL "0")
            message(FATAL_ERROR "pybind11 not installed, cannot compile Python bindings!")
        endif()
        find_package(pybind11 REQUIRED HINTS "${_PYBIND11_CMAKEDIR}")
    endif()

    pybind11_add_module(
        ${PYLIB}
        MODULE
        NO_EXTRAS  # No -flto, it makes linking stage annoyingly slow.
        "${SRC_DIR}/operators/ImportExportUniform.cpp"
        "${SRC_DIR}/bindings/Bindings.cpp"
        "${SRC_DIR}/bindings/Fields.cpp"
        "${SRC_DIR}/bindings/Obstacles.cpp"
        "${SRC_DIR}/bindings/Operators.cpp"
        "${SRC_DIR}/bindings/Simulation.cpp"
    )
    # Manually add pybind11 include dirs before those from Python.
    target_include_directories(${PYLIB} BEFORE PUBLIC ${pybind11_INCLUDE_DIRS})
    target_link_libraries(${PYLIB} PUBLIC ${CORE})
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
