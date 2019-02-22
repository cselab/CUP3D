//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#include "PressurePenalization.h"

#ifdef _ACCFFT_
#include "../poisson/PoissonSolverACCPeriodic.h"
#include "../poisson/PoissonSolverACCUnbounded.h"
#else
#include "../poisson/PoissonSolverPeriodic.h"
#include "../poisson/PoissonSolverUnbounded.h"
#endif
// TODO : Cosine transform on GPU!?
#include "../poisson/PoissonSolverMixed.h"
#include "../poisson/PoissonSolverHYPREMixed.h"

class KernelGradP
{
 private:
  const double dt;
  const Real extent[3];

 public:
  const std::array<int, 3> stencil_start = {-1, -1, -1};
  const std::array<int, 3> stencil_end = {2, 2, 2};
  const StencilInfo stencil = StencilInfo(-1,-1,-1, 2,2,2, false, 1, 4);

  KernelGradP(double _dt,const Real ext[3]): dt(_dt), extent{ext[0],ext[1],ext[2]} {}

  ~KernelGradP() {}

  template <typename Lab, typename BlockType>
  void operator()(Lab & lab, const BlockInfo& info, BlockType& o) const
  {
    const Real fac = - 0.5 * dt / info.h_gridpoint;
    for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
    for(int iy=0; iy<FluidBlock::sizeY; ++iy)
    for(int ix=0; ix<FluidBlock::sizeX; ++ix)
    {
      // p contains the pressure correction after the Poisson solver
      const Real Uf = o(ix,iy,iz).u + fac*(lab(ix+1,iy,iz).p-lab(ix-1,iy,iz).p);
      const Real Vf = o(ix,iy,iz).v + fac*(lab(ix,iy+1,iz).p-lab(ix,iy-1,iz).p);
      const Real Wf = o(ix,iy,iz).w + fac*(lab(ix,iy,iz+1).p-lab(ix,iy,iz-1).p);
      const Real US=o(ix,iy,iz).tmpU, VS=o(ix,iy,iz).tmpV, WS=o(ix,iy,iz).tmpW;
      o(ix,iy,iz).u = ( Uf + o(ix,iy,iz).chi * US) / (1+o(ix,iy,iz).chi);
      o(ix,iy,iz).v = ( Vf + o(ix,iy,iz).chi * VS) / (1+o(ix,iy,iz).chi);
      o(ix,iy,iz).w = ( Wf + o(ix,iy,iz).chi * WS) / (1+o(ix,iy,iz).chi);
    }
  }
};

PressurePenalization::PressurePenalization(SimulationData & s) : Operator(s)
{
  if(sim.bUseFourierBC)
  pressureSolver = new PoissonSolverPeriodic(sim);
  else if (sim.bUseUnboundedBC)
  pressureSolver = new PoissonSolverUnbounded(sim);
  else pressureSolver = new PoissonSolverMixed(sim);
  sim.pressureSolver = pressureSolver;
}

void PressurePenalization::operator()(const double dt)
{
  pressureSolver->solve();

  sim.startProfiler("GradP Penal");
  { //pressure correction dudt* = - grad P / rho
    const int nthreads = omp_get_max_threads();
    std::vector<KernelGradP*> diff(nthreads, nullptr);
    #pragma omp parallel for schedule(static, 1)
    for(int i=0;i<nthreads;++i) diff[i] = new KernelGradP(dt, sim.extent);

    compute<KernelGradP>(diff);
    for(int i=0; i<nthreads; i++) delete diff[i];
  }
  sim.stopProfiler();

  check("pressure - end");
}
