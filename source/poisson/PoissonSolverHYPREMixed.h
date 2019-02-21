//
//  CubismUP_2D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//


#pragma once
#ifdef CUP_HYPRE

#include "HYPRE_struct_ls.h"
#include "PoissonSolver.h"

class PoissonSolverMixed_HYPRE : public PoissonSolver
{
  const std::string solver;
  HYPRE_StructGrid     hypre_grid;
  HYPRE_StructStencil  hypre_stencil;
  HYPRE_StructMatrix   hypre_mat;
  HYPRE_StructVector   hypre_rhs;
  HYPRE_StructVector   hypre_sol;
  HYPRE_StructSolver   hypre_solver;
  HYPRE_StructSolver   hypre_precond;
  Real pLast = 0;

  int peidx_0() const {  int ret[3]; grid.peindex(ret); return ret[0]; }
  int peidx_1() const {  int ret[3]; grid.peindex(ret); return ret[1]; }
  int peidx_2() const {  int ret[3]; grid.peindex(ret); return ret[2]; }
  const int peidx[3] = {peidx_0(), peidx_1(), peidx_2()};
  HYPRE_Int ilower[] = {
    (int) myN[0] * peidx[0],
    (int) myN[1] * peidx[1],
    (int) myN[2] * peidx[2]
  };
  HYPRE_Int iupper[] = {
    (int) myN[0] * (peidx[0]+1) - 1,
    (int) myN[1] * (peidx[1]+1) - 1,
    (int) myN[2] * (peidx[2]+1) - 1
  };

  inline size_t linaccess(const size_t ix,const size_t iy,const size_t iz) const
  {
    return ix + myN[0] * iy + myN[0]*myN[1] * iz;
  }

 public:
  void solve() override;

  PoissonSolverMixed_HYPRE(SimulationData& s);

  std::string getName() {
    return "hypre";
  }

  ~PoissonSolverMixed_HYPRE();
};

#endif
