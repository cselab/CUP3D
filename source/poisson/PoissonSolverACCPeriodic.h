//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#ifndef CubismUP_3D_PoissonSolverPeriodicACC_h
#define CubismUP_3D_PoissonSolverPeriodicACC_h

#include "poisson/PoissonSolver.h"

CubismUP_3D_NAMESPACE_BEGIN

class PoissonSolverPeriodic : public PoissonSolver
{
  MPI_Comm c_comm;
  // the local pencil size and the allocation size
  int isize[3], osize[3], istart[3], ostart[3];
  const size_t gz_hat = gsize[2] / 2 + 1;
  size_t alloc_max;
  //Real * rho_gpu;
  //Real * phi_gpu;
  Real * phi_hat;
  void * plan;

  const double h = sim.uniformH();

public:
  PoissonSolverPeriodic(SimulationData & s);

  void solve();
  void testComm();
  ~PoissonSolverPeriodic();
};

CubismUP_3D_NAMESPACE_END
#endif // CubismUP_3D_PoissonSolverPeriodicACC_h
