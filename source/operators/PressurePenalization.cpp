//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#include "PressurePenalization.h"
#include "PenalizationObstacleVisitor.h"

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
#include "../poisson/PoissonSolverPETSCMixed.h"

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
      const Real Uf = o(ix,iy,iz).u + fac*(lab(ix+1,iy,iz).p-lab(ix-1,iy,iz).p);
      const Real Vf = o(ix,iy,iz).v + fac*(lab(ix,iy+1,iz).p-lab(ix,iy-1,iz).p);
      const Real Wf = o(ix,iy,iz).w + fac*(lab(ix,iy,iz+1).p-lab(ix,iy,iz-1).p);
      const Real US=o(ix,iy,iz).tmpU, VS=o(ix,iy,iz).tmpV, WS=o(ix,iy,iz).tmpW;
      #if PENAL_TYPE==0
       o(ix,iy,iz).u = Uf + o(ix,iy,iz).chi * US; // explicit penal part 2
       o(ix,iy,iz).v = Vf + o(ix,iy,iz).chi * VS; // explicit penal part 2
       o(ix,iy,iz).w = Wf + o(ix,iy,iz).chi * WS; // (part 1 in advectdiffuse)
      #else // implicit
       // p contains the pressure correction after the Poisson solver
       o(ix,iy,iz).u = ( Uf + o(ix,iy,iz).chi * US) / (1+o(ix,iy,iz).chi);
       o(ix,iy,iz).v = ( Vf + o(ix,iy,iz).chi * VS) / (1+o(ix,iy,iz).chi);
       o(ix,iy,iz).w = ( Wf + o(ix,iy,iz).chi * WS) / (1+o(ix,iy,iz).chi);
      #endif
    }
  }
};

/*
class KernelGradP_nonUniform
{
  const Real dt, extent[3];
 public:
  const std::array<int, 3> stencil_start = {-1, -1, -1};
  const std::array<int, 3> stencil_end = {2, 2, 2};
  const StencilInfo stencil = StencilInfo(-1,-1,-1, 2,2,2, false, 1, 4);

  KernelGradP_nonUniform(double _dt,const Real ext[3]): dt(_dt), extent{ext[0],ext[1],ext[2]} {}

  template <typename Lab, typename BlockType>
  void operator()(Lab & lab, const BlockInfo& info, BlockType& o) const
  {
    // FD coefficients for first derivative
    const BlkCoeffX& cx = o.fd_cx.first;
    const BlkCoeffY& cy = o.fd_cy.first;
    const BlkCoeffZ& cz = o.fd_cz.first;
    for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
    for(int iy=0; iy<FluidBlock::sizeY; ++iy)
    for(int ix=0; ix<FluidBlock::sizeX; ++ix)
    {
      const FluidElement &L =lab(ix,iy,iz);
      const FluidElement &LW=lab(ix-1,iy,iz), &LE=lab(ix+1,iy,iz);
      const FluidElement &LS=lab(ix,iy-1,iz), &LN=lab(ix,iy+1,iz);
      const FluidElement &LF=lab(ix,iy,iz-1), &LB=lab(ix,iy,iz+1);
      const Real Uf = o(ix,iy,iz).u - dt * __FD_2ND(ix, cx, LW.p, L.p, LE.p);
      const Real Vf = o(ix,iy,iz).v - dt * __FD_2ND(iy, cy, LS.p, L.p, LN.p);
      const Real Wf = o(ix,iy,iz).w - dt * __FD_2ND(iz, cz, LF.p, L.p, LB.p);
      const Real US=o(ix,iy,iz).tmpU, VS=o(ix,iy,iz).tmpV, WS=o(ix,iy,iz).tmpW;
      #if PENAL_TYPE==0
       o(ix,iy,iz).u = Uf + o(ix,iy,iz).chi * US; // explicit penal part 2
       o(ix,iy,iz).v = Vf + o(ix,iy,iz).chi * VS; // explicit penal part 2
       o(ix,iy,iz).w = Wf + o(ix,iy,iz).chi * WS; // (part 1 in advectdiffuse)
      #else // implicit
       // p contains the pressure correction after the Poisson solver
       o(ix,iy,iz).u = ( Uf + o(ix,iy,iz).chi * US) / (1+o(ix,iy,iz).chi);
       o(ix,iy,iz).v = ( Vf + o(ix,iy,iz).chi * VS) / (1+o(ix,iy,iz).chi);
       o(ix,iy,iz).w = ( Wf + o(ix,iy,iz).chi * WS) / (1+o(ix,iy,iz).chi);
      #endif
    }
  }
};
*/

PressurePenalization::PressurePenalization(SimulationData & s) : Operator(s)
{
  if(sim.bUseFourierBC)
  pressureSolver = new PoissonSolverPeriodic(sim);
  else if (sim.bUseUnboundedBC)
  pressureSolver = new PoissonSolverUnbounded(sim);
  #ifdef CUP_HYPRE
  else if (sim.useSolver == "hypre")
  pressureSolver = new PoissonSolverMixed_HYPRE(sim);
  #endif
  #ifdef CUP_PETSC
  else if (sim.useSolver == "petsc")
  pressureSolver = new PoissonSolverMixed_PETSC(sim);
  #endif
  else
  pressureSolver = new PoissonSolverMixed(sim);
  sim.pressureSolver = pressureSolver;
}

void PressurePenalization::operator()(const double dt)
{
  pressureSolver->solve();

  #if PENAL_TYPE==0
  sim.startProfiler("PresRHS Uobst.");
  {
    //zero fields, going to contain Udef:
    #pragma omp parallel for schedule(static)
    for(unsigned i=0; i<vInfo.size(); i++) {
      FluidBlock& b = *(FluidBlock*)vInfo[i].ptrBlock;
      for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
      for(int iy=0; iy<FluidBlock::sizeY; ++iy)
      for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
        b(ix,iy,iz).tmpU = 0; b(ix,iy,iz).tmpV = 0; b(ix,iy,iz).tmpW = 0;
      }
    }
    //store deformation velocities onto tmp fields:
    ObstacleVisitor*visitor=new PenalizationObstacleVisitor(grid,dt,sim.uinf);
    sim.obstacle_vector->Accept(visitor);
    delete visitor;
  }
  sim.stopProfiler();
  #endif

  sim.startProfiler("GradP Penal");
  { //pressure correction dudt* = - grad P / rho
    const KernelGradP K(dt, sim.extent);
    compute<KernelGradP>(K);
  }
  sim.stopProfiler();

  check("pressure - end");
}
