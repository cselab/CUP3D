//
//  Cubism3D
//  Copyright (c) 2019 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Ivica Kicic (kicici@ethz.ch).
//

#pragma once

/*
 * Macros and definition used by other header files.
 */

// Are we compiling from CubismUP_3D's makefile?
#ifndef CUP_NO_MACROS_HEADER
// No, it's either CMake or external code. Load compile-time settings from this header file.
#include <generated.CubismUP3DMacros.h>
#endif

#ifndef CubismUP_3D_NAMESPACE_BEGIN
#define CubismUP_3D_NAMESPACE_BEGIN namespace cubismup3d {
#endif

#ifndef CubismUP_3D_NAMESPACE_END
#define CubismUP_3D_NAMESPACE_END   }  // namespace cubismup3d
#endif

#ifndef CUP_ALIGNMENT
#define CUP_ALIGNMENT 64
#endif
#define CUBISM_ALIGNMENT CUP_ALIGNMENT

#ifdef _FLOAT_PRECISION_
using Real = float;
#define MPI_Real MPI_FLOAT
#endif
#ifdef _DOUBLE_PRECISION_
using Real = double;
#define MPI_Real MPI_DOUBLE
#endif
#ifdef _LONG_DOUBLE_PRECISION_
using Real = long double;
#define MPI_Real MPI_LONG_DOUBLE
#endif
