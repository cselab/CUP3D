//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Hugues de Laroussilhe.
//

#include "SpectralManipACC.h"
#include <cuda_runtime_api.h>
#include "../poisson/PoissonSolverACC_common.h"
#include "accfft_common.h"
#ifndef CUP_SINGLE_PRECISION
  #include "accfft_gpu.h"
  typedef accfft_plan_gpu acc_plan;
#else
  #include "accfft_gpuf.h"
  typedef accfft_plan_gpuf acc_plan;
#endif

#ifndef CUP_SINGLE_PRECISION
#define MPIREAL MPI_DOUBLE
#else
#define MPIREAL MPI_FLOAT
#endif /* CUP_SINGLE_PRECISION */

CubismUP_3D_NAMESPACE_BEGIN
using namespace cubism;

void SpectralManipPeriodic::_compute_largeModesForcing()
{
}

void SpectralManipPeriodic::_compute_analysis()
{
}

void SpectralManipPeriodic::_compute_IC(const std::vector<Real> &K,
                                        const std::vector<Real> &E)
{
}

SpectralManipPeriodic::SpectralManipPeriodic(SimulationData&s): SpectralManip(s)
{
}

void SpectralManipPeriodic::prepareFwd()
{
}

void SpectralManipPeriodic::prepareBwd()
{
}

void SpectralManipPeriodic::runFwd() const
{
}

void SpectralManipPeriodic::runBwd() const
{
}

SpectralManipPeriodic::~SpectralManipPeriodic()
{
}

CubismUP_3D_NAMESPACE_END
#undef MPIREAL
