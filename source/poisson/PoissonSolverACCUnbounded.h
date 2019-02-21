//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#pragma once

#include "PoissonSolver.h"

class PoissonSolverUnbounded : public PoissonSolver
{
  // the local pencil size and the allocation size
  int isize[3], osize[3], istart[3], ostart[3];
  const int mx = 2*gsize[0]-1, my = 2*gsize[1]-1, mz = 2*gsize[2]-1;
  const size_t mz_pad = mz/2 +1, myftNx = (mx+1)/m_size;
  const int szFft[3] = {(int) myftNx, (int) gsize[1], (int) gsize[2] };
  const int szCup[3] = {std::min(szFft[0],(int)myN[0]),(int)myN[1],(int)myN[2]};

  MPI_Comm c_comm;
  size_t alloc_max;
  Real * fft_rhs;
  Real * gpuGhat;
  Real * gpu_rhs;
  void * plan;
  MPI_Datatype submat;

public:
  PoissonSolverUnbounded(SimulationData & s);

  void solve();

  void cub2padded() const;
  void padded2cub() const;
  void padded2gpu() const;
  void gpu2padded() const;
  ~PoissonSolverUnbounded();
};
