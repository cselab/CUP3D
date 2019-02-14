//
//  CubismUP_2D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#pragma once

#include "../SimulationData.h"
#include <fftw3-mpi.h>

#include "Cubism/BlockInfo.h"

#include <vector>
#include <cassert>
#include <cstring>

#ifndef CUP_SINGLE_PRECISION
#define _FFTW_(s) fftw_##s
typedef fftw_complex mycomplex;
typedef fftw_plan myplan;
#define MPIREAL MPI_DOUBLE
#else
#define _FFTW_(s) fftwf_##s
typedef fftwf_complex mycomplex;
typedef fftwf_plan myplan;
#define MPIREAL MPI_FLOAT
#endif /* CUP_SINGLE_PRECISION */

class PoissonSolver
{
 public:
    PoissonSolver() {}
    PoissonSolver(const PoissonSolver& c) = delete;
    virtual ~PoissonSolver() {}

    virtual void solve() = 0;
};
