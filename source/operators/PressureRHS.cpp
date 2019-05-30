//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#include "operators/PressureRHS.h"
#include "poisson/PoissonSolver.h"
#include "obstacles/ObstacleVector.h"

CubismUP_3D_NAMESPACE_BEGIN
using namespace cubism;

namespace {

using CHIMAT = Real[CUP_BLOCK_SIZE][CUP_BLOCK_SIZE][CUP_BLOCK_SIZE];
using UDEFMAT = Real[CUP_BLOCK_SIZE][CUP_BLOCK_SIZE][CUP_BLOCK_SIZE][3];

template<int withObstacles = 1>
class KernelPressureRHS
{
  const Real dt;
  PoissonSolver * const solver;

 public:
  const std::array<int, 3> stencil_start = {-1,-1,-1}, stencil_end = {2, 2, 2};
  const StencilInfo stencil =
    withObstacles ? StencilInfo(-1,-1,-1, 2,2,2, false, 6, 1,2,3, 5,6,7)
                  : StencilInfo(-1,-1,-1, 2,2,2, false, 6, 1,2,3);

  KernelPressureRHS(double _dt, PoissonSolver* ps) : dt(_dt), solver(ps) {}

  template <typename Lab, typename BlockType>
  void operator()(Lab & lab, const BlockInfo& info, BlockType& o) const
  {
    const Real h = info.h_gridpoint, fac = .5*h*h/dt;
    Real* __restrict__ const ret = solver->data + solver->_offset_ext(info);
    const unsigned SX=solver->stridex, SY=solver->stridey, SZ=solver->stridez;

    for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
    for(int iy=0; iy<FluidBlock::sizeY; ++iy)
    for(int ix=0; ix<FluidBlock::sizeX; ++ix)
    {
      const FluidElement &L  = lab(ix,  iy,  iz);
      const FluidElement &LW = lab(ix-1,iy,  iz  ), &LE = lab(ix+1,iy,  iz  );
      const FluidElement &LS = lab(ix,  iy-1,iz  ), &LN = lab(ix,  iy+1,iz  );
      const FluidElement &LF = lab(ix,  iy,  iz-1), &LB = lab(ix,  iy,  iz+1);
      if(withObstacles) {
        const Real divUs = LE.tmpU-LW.tmpU + LN.tmpV-LS.tmpV + LB.tmpW-LF.tmpW;
        const Real divUf = LE.u-LW.u + LN.v-LS.v + LB.w-LF.w;
        ret[SZ*iz +SY*iy +SX*ix] = fac*(divUf - L.chi*divUs);
      } else
        ret[SZ*iz +SY*iy +SX*ix] = fac*(LE.u-LW.u + LN.v-LS.v + LB.w-LF.w);
      //o(ix,iy,iz).p = ret[SZ*iz + SY*iy + SX*ix];
    }
  }
};

template<int withObstacles = 1>
class KernelPressureRHS_nonUniform
{
 private:
  const Real dt, invdt;
  PoissonSolver * const solver;

 public:
  const std::array<int, 3> stencil_start = {-1,-1,-1}, stencil_end = {2, 2, 2};
  const StencilInfo stencil =
    withObstacles ? StencilInfo(-1,-1,-1, 2,2,2, false, 6, 1,2,3, 5,6,7)
                  : StencilInfo(-1,-1,-1, 2,2,2, false, 6, 1,2,3);

  KernelPressureRHS_nonUniform(double _dt, PoissonSolver* ps) :
    dt(_dt), invdt(1/_dt), solver(ps) {}

  template <typename Lab, typename BlockType>
  void operator()(Lab & lab, const BlockInfo& info, BlockType& o) const
  {
    // FD coefficients for first derivative
    const BlkCoeffX &cx =o.fd_cx.first, &cy =o.fd_cy.first, &cz =o.fd_cz.first;
    Real* __restrict__ const ret = solver->data + solver->_offset_ext(info);
    const unsigned SX=solver->stridex, SY=solver->stridey, SZ=solver->stridez;

    for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
    for(int iy=0; iy<FluidBlock::sizeY; ++iy)
    for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
      Real h[3]; info.spacing(h, ix, iy, iz);
      const Real fac = h[0]*h[1]*h[2]*invdt;
      const FluidElement& L  = lab(ix,  iy,  iz);
      const FluidElement& LW = lab(ix-1,iy,  iz  ), & LE = lab(ix+1,iy,  iz  );
      const FluidElement& LS = lab(ix,  iy-1,iz  ), & LN = lab(ix,  iy+1,iz  );
      const FluidElement& LF = lab(ix,  iy,  iz-1), & LB = lab(ix,  iy,  iz+1);
      const Real dudx = __FD_2ND(ix, cx, LW.u, L.u, LE.u);
      const Real dvdy = __FD_2ND(iy, cy, LS.v, L.v, LN.v);
      const Real dwdz = __FD_2ND(iz, cz, LF.w, L.w, LB.w);
      if(withObstacles) {
        const Real dusdx = __FD_2ND(ix, cx, LW.tmpU, L.tmpU, LE.tmpU);
        const Real dvsdy = __FD_2ND(iy, cy, LS.tmpV, L.tmpV, LN.tmpV);
        const Real dwsdz = __FD_2ND(iz, cz, LF.tmpW, L.tmpW, LB.tmpW);
        const Real divUs = dusdx + dvsdy + dwsdz, divUf = dudx + dvdy + dwdz;
        ret[SZ*iz + SY*iy + SX*ix] = fac * (divUf - L.chi * divUs);
      } else ret[SZ*iz + SY*iy + SX*ix] = fac * (dudx + dvdy + dwdz);
      //o(ix,iy,iz).p = fac * RHS(lab, ix,iy,iz, pFac);
    }
  }
};

struct PressureRHSObstacleVisitor : public ObstacleVisitor
{
  FluidGridMPI * const grid;
  const std::vector<cubism::BlockInfo>& vInfo = grid->getBlocksInfo();

  PressureRHSObstacleVisitor(FluidGridMPI*g) : grid(g) { }

  void visit(Obstacle* const obstacle)
  {
    #pragma omp parallel
    {
      const auto& obstblocks = obstacle->getObstacleBlocks();
      #pragma omp for schedule(dynamic, 1)
      for (size_t i = 0; i < vInfo.size(); ++i)
      {
        const cubism::BlockInfo& info = vInfo[i];
        const auto pos = obstblocks[info.blockID];
        if(pos == nullptr) continue;

        FluidBlock& b = *(FluidBlock*)info.ptrBlock;
        const UDEFMAT & __restrict__ UDEF = pos->udef;
        const CHIMAT & __restrict__ CHI = pos->chi;
        //const CHIMAT & __restrict__ SDF = pos->sdf;

        for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
        for(int iy=0; iy<FluidBlock::sizeY; ++iy)
        for(int ix=0; ix<FluidBlock::sizeX; ++ix)
        {
          // What if multiple obstacles share a block? Do not write udef onto
          // grid if CHI stored on the grid is greater than obst's CHI.
          if(b(ix,iy,iz).chi > CHI[iz][iy][ix]) continue;
          // What if two obstacles overlap? Let's plus equal. After all here
          // we are computing divUs, maybe one obstacle has divUs 0. We will
          // need a repulsion term of the velocity at some point in the code.
          b(ix,iy,iz).tmpU += UDEF[iz][iy][ix][0];
          b(ix,iy,iz).tmpV += UDEF[iz][iy][ix][1];
          b(ix,iy,iz).tmpW += UDEF[iz][iy][ix][2];
        }
      }
    }
  }
};

}

void PressureRHS::operator()(const double dt)
{
  sim.startProfiler("PresRHS Udef");
  if(sim.obstacle_vector->nObstacles() > 0)
  { //zero fields, going to contain Udef:
    #pragma omp parallel for schedule(static)
    for(size_t i=0; i<vInfo.size(); i++) {
      FluidBlock& b = *(FluidBlock*)vInfo[i].ptrBlock;
      for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
      for(int iy=0; iy<FluidBlock::sizeY; ++iy)
      for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
        b(ix,iy,iz).tmpU = 0; b(ix,iy,iz).tmpV = 0; b(ix,iy,iz).tmpW = 0;
      }
    }
    //store deformation velocities onto tmp fields:
    ObstacleVisitor* visitor = new PressureRHSObstacleVisitor(grid);
    sim.obstacle_vector->Accept(visitor);
    delete visitor;
  }
  sim.stopProfiler();

  sim.startProfiler("PresRHS Kernel");
  //place onto p: ( div u^(t+1) - div u^* ) / dt
  //where i want div u^(t+1) to be equal to div udef
  sim.pressureSolver->reset();

  if(sim.bUseStretchedGrid)
  {
    if(sim.obstacle_vector->nObstacles())
    {
      const KernelPressureRHS_nonUniform<1> K(dt, sim.pressureSolver);
      compute<KernelPressureRHS_nonUniform<1>>(K);
    }
    else
    {
      const KernelPressureRHS_nonUniform<0> K(dt, sim.pressureSolver);
      compute<KernelPressureRHS_nonUniform<0>>(K);
    }
  }
  else
  {
    if(sim.obstacle_vector->nObstacles())
    {
      const KernelPressureRHS<1> K(dt, sim.pressureSolver);
      compute<KernelPressureRHS<1>>(K);
    }
    else
    {
      const KernelPressureRHS<0> K(dt, sim.pressureSolver);
      compute<KernelPressureRHS<0>>(K);
    }
  }
  sim.stopProfiler();

  check("PressureRHS");
}

CubismUP_3D_NAMESPACE_END
