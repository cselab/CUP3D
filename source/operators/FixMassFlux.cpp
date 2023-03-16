//
//  Cubism3D
//  Copyright (c) 2023 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//

#include "FixMassFlux.h"

CubismUP_3D_NAMESPACE_BEGIN using namespace cubism;

static Real avgUx_nonUniform(const std::vector<BlockInfo>& myInfo, const Real* const uInf, const Real volume)
{
  // Average Ux on the simulation volume :
  //   Sum on the xz-plane (uniform)
  //   Integral along Y    (non-uniform)
  //
  // <Ux>_{xz} (iy) = 1/(Nx.Ny) . \Sum_{ix, iz} u(ix,iy,iz)
  //
  // <Ux>_{xyz} = 1/Ly . \Sum_{iy} <Ux>_{xz} (iy) . h_y(*,iy,*)
  //            = /1(Nx.Ny.Ly) . \Sum_{ix,iy,iz} u(ix,iy,iz).h_y(ix,iy,iz)
  Real avgUx = 0.;
  const int nBlocks = myInfo.size();

  #pragma omp parallel for reduction(+ : avgUx)
  for (int i = 0; i < nBlocks; i++)
  {
    const BlockInfo& info = myInfo[i];
    const VectorBlock& b = *(const VectorBlock*)info.ptrBlock;
    const Real h3 = info.h*info.h*info.h;
    for (int z = 0; z < VectorBlock::sizeZ; ++z)
    for (int y = 0; y < VectorBlock::sizeY; ++y)
    for (int x = 0; x < VectorBlock::sizeX; ++x)
    {
      avgUx += (b(x,y,z).u[0] + uInf[0]) * h3;
    }
  }
  avgUx = avgUx / volume;
  return avgUx;
}

FixMassFlux::FixMassFlux(SimulationData& s): Operator(s) {}

void FixMassFlux::operator()(const double dt)
{
  sim.startProfiler("FixedMassFlux");

  const std::vector<BlockInfo> &  velInfo = sim.velInfo();

  // fix base_u_avg and y_max AD HOC for channel flow
  const Real volume = sim.extents[0]*sim.extents[1]*sim.extents[2];
  const Real y_max = sim.extents[1];
  const Real u_avg = 2.0/3.0 * sim.uMax_forced;
  Real u_avg_msr = avgUx_nonUniform(velInfo, sim.uinf.data(), volume);
  MPI_Allreduce(MPI_IN_PLACE, &u_avg_msr, 1, MPI_Real, MPI_SUM, sim.comm);
  const Real delta_u = u_avg - u_avg_msr;
  const Real reTau = std::sqrt(std::fabs(delta_u/sim.dt)) / sim.nu;
  const Real scale = 6*delta_u;

  if (sim.rank == 0)
  {
    printf(
        "Measured <Ux>_V = %25.16e,\n"
        "target   <Ux>_V = %25.16e,\n"
        "delta    <Ux>_V = %25.16e,\n"
        "scale           = %25.16e,\n"
        "Re_tau          = %25.16e,\n",
        u_avg_msr, u_avg, delta_u, scale, reTau);
  }
 
  #pragma omp parallel for
  for(size_t i=0; i<velInfo.size(); i++)
  {
      VectorBlock& v = *(VectorBlock*)velInfo[i].ptrBlock;
      for (int z = 0; z < VectorBlock::sizeZ; ++z)
      for (int y = 0; y < VectorBlock::sizeY; ++y)
      {
           Real p[3];
           velInfo[i].pos(p, 0, y, 0);
           const Real aux = 6 * scale * p[1]/y_max * (1.0 - p[1]/y_max);
           for (int x = 0; x < VectorBlock::sizeX; ++x)
                v(x,y,z).u[0] += aux;
      }
  }

  sim.stopProfiler();
}

CubismUP_3D_NAMESPACE_END
