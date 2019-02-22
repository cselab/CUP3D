//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#pragma once

#include "PoissonSolver.h"

class PoissonSolverPeriodic : public PoissonSolver
{
  void * fwd, * bwd;
  const size_t nz_hat = gsize[2]/2+1;
  const double norm_factor = 1./(gsize[0]*h*gsize[1]*h*gsize[2]*h);
  ptrdiff_t alloc_local=0, local_n0=0, local_0_start=0, local_n1=0, local_1_start=0;

 protected:

  void _solve();

 public:

  PoissonSolverPeriodic(SimulationData & s);

  void solve();

  ~PoissonSolverPeriodic();
};
