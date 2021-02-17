//
//  CubismUP_3D
//  Copyright (c) 2020 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Michalis Chatzimanolakis (michaich@ethz.ch).
//

#include "PressureProjection.h"

CubismUP_3D_NAMESPACE_BEGIN
using namespace cubism;

namespace {

class KernelGradP
{
  const Real dt;
 public:
  const std::array<int, 3> stencil_start = {-1,-1,-1}, stencil_end = {2, 2, 2};
  const StencilInfo stencil{-1,-1,-1, 2,2,2, false, {{FE_P}}};

  KernelGradP(double _dt, const std::array<Real, 3> &ext): dt(_dt) {}

  ~KernelGradP() {}

  template <typename Lab, typename BlockType>
  void operator()(Lab & lab, const BlockInfo& info, BlockType& o) const
  {
    const Real fac = - 0.5 * dt / info.h_gridpoint;
    for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
    for(int iy=0; iy<FluidBlock::sizeY; ++iy)
    for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
       // p contains the pressure correction after the Poisson solver
       o(ix,iy,iz).u += fac*(lab(ix+1,iy,iz).p-lab(ix-1,iy,iz).p);
       o(ix,iy,iz).v += fac*(lab(ix,iy+1,iz).p-lab(ix,iy-1,iz).p);
       o(ix,iy,iz).w += fac*(lab(ix,iy,iz+1).p-lab(ix,iy,iz-1).p);
    }
  }
};

}

PressureProjection::PressureProjection(SimulationData & s) : Operator(s)
{
  pressureSolver = new PoissonSolverAMR(sim);
  sim.pressureSolver = pressureSolver;
}

void PressureProjection::operator()(const double dt)
{
  const std::vector<cubism::BlockInfo>& vInfo = sim.vInfo();

  if (sim.TimeOrder == 2 && sim.step >= sim.step_2nd_start-1)
  {
    #pragma omp parallel for
    for(size_t i=0; i<vInfo.size(); i++)
    {
      FluidBlock& b = *(FluidBlock*)vInfo[i].ptrBlock;
      for (int iz=0; iz<FluidBlock::sizeZ; ++iz)
      for (int iy=0; iy<FluidBlock::sizeY; ++iy)
      for (int ix=0; ix<FluidBlock::sizeX; ++ix)
      {
        b.dataOld[iz][iy][ix][3] = b(ix,iy,iz).p;
        b(ix,iy,iz).p = 0.0;
      }
    }
  }

  pressureSolver->solve();

  sim.startProfiler("GradP"); //pressure correction dudt* = - grad P / rho
  const KernelGradP K( (sim.TimeOrder == 1 || sim.step < sim.step_2nd_start) ? dt:(dt/sim.coefU[0]), sim.extent);
  compute<KernelGradP>(K);
  sim.stopProfiler();

  if (sim.TimeOrder == 2 && sim.step >= sim.step_2nd_start)
  {
    #pragma omp parallel for
    for(size_t i=0; i<vInfo.size(); i++)
    {
      FluidBlock& b = *(FluidBlock*)vInfo[i].ptrBlock;
      for (int iz=0; iz<FluidBlock::sizeZ; ++iz)
      for (int iy=0; iy<FluidBlock::sizeY; ++iy)
      for (int ix=0; ix<FluidBlock::sizeX; ++ix)
        b(ix,iy,iz).p += b.dataOld[iz][iy][ix][3];
    }
  }

  check("PressureProjection");
}

CubismUP_3D_NAMESPACE_END
