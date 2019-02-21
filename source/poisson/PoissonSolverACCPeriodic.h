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
  MPI_Comm c_comm;
  // the local pencil size and the allocation size
  int isize[3], osize[3], istart[3], ostart[3];
  size_t alloc_max;
  Real * rho_gpu;
  Real * phi_gpu;
  Real * phi_hat;
  void * plan;

public:
  PoissonSolverPeriodic(SimulationData & s);

  void solve();

  ~PoissonSolverPeriodic();
};
