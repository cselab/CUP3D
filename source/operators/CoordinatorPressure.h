//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch) and Christian Conti.
//

#ifndef CubismUP_3D_CoordinatorPressure_h
#define CubismUP_3D_CoordinatorPressure_h

#include "PenalizationObstacleVisitor.h"
#include "../obstacles/IF3D_ObstacleVector.h"

//#include "../poisson/PoissonSolverScalarACC_freespace.h"
//#include "../poisson/PoissonSolverScalarACC.h"

#include "../poisson/PoissonSolverUnbounded.h"
#include "../poisson/PoissonSolverPeriodic.h"
#include "../poisson/PoissonSolverMixed.h"

class OperatorDivergence : public GenericLabOperator
{
 private:
  double dt;
  const Real ext[3], fadeLen[3];
  const Real iFade[3] = {1/fadeLen[0], 1/fadeLen[1], 1/fadeLen[2]};

  inline bool _is_touching(const BlockInfo& i) const
  {
    Real maxP[3], minP[3]; i.pos(minP, 0, 0, 0); const Real& h = i.h_gridpoint;
    i.pos(maxP, CUP_BLOCK_SIZE-1, CUP_BLOCK_SIZE-1, CUP_BLOCK_SIZE-1);
    const bool touchW=h+fadeLen[0]>=minP[0],touchE=h+fadeLen[0]>=ext[0]-maxP[0];
    const bool touchS=h+fadeLen[1]>=minP[1],touchN=h+fadeLen[1]>=ext[1]-maxP[1];
    const bool touchB=h+fadeLen[2]>=minP[2],touchF=h+fadeLen[2]>=ext[2]-maxP[2];
    return touchN || touchE || touchS || touchW || touchF || touchB;
  }

  inline Real fade(const BlockInfo&i, const int x,const int y,const int z) const
  {
    Real p[3]; i.pos(p, x, y, z); const Real& h = i.h_gridpoint;
    const Real zt = iFade[2] * std::max(Real(0), h+fadeLen[2] -(ext[2]-p[2]) );
    const Real zb = iFade[2] * std::max(Real(0), h+fadeLen[2] - p[2] );
    const Real yt = iFade[1] * std::max(Real(0), h+fadeLen[1] -(ext[1]-p[1]) );
    const Real yb = iFade[1] * std::max(Real(0), h+fadeLen[1] - p[1] );
    const Real xt = iFade[0] * std::max(Real(0), h+fadeLen[0] -(ext[0]-p[0]) );
    const Real xb = iFade[0] * std::max(Real(0), h+fadeLen[0] - p[0] );
    return 1-std::pow(std::min( std::max({zt,zb,yt,yb,xt,xb}), (Real)1), 2);
  }

 public:
  OperatorDivergence(double _dt, const Real buf[3], const Real extent[3])
   : dt(_dt), ext{extent[0],extent[1],extent[2]}, fadeLen{buf[0],buf[1],buf[2]}
  {
    stencil_start[0] = -1; stencil_start[1] = -1;  stencil_start[2] = -1;
    stencil_end[0] = 2;  stencil_end[1] = 2;  stencil_end[2] = 2;
    stencil = StencilInfo(-1,-1,-1, 2,2,2, false, 3, 0,1,2,3);
  }
  ~OperatorDivergence() {}

  template <typename Lab, typename BlockType>
  void operator()(Lab & lab, const BlockInfo& info, BlockType& o) const
  {
    const Real h = info.h_gridpoint, fac = .5*h*h/dt;
    if( not _is_touching(info) )
    {
      for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
      for(int iy=0; iy<FluidBlock::sizeY; ++iy)
      for(int ix=0; ix<FluidBlock::sizeX; ++ix)
      {
        const FluidElement &L =lab(ix,iy,iz);
        const FluidElement &LW=lab(ix-1,iy,iz), &LE=lab(ix+1,iy,iz);
        const FluidElement &LS=lab(ix,iy-1,iz), &LN=lab(ix,iy+1,iz);
        const FluidElement &LF=lab(ix,iy,iz-1), &LB=lab(ix,iy,iz+1);
        const Real dU = LE.u-LW.u, dV = LN.v-LS.v, dW = LB.w-LF.w;
        const Real dXx_ux = (LE.chi-LW.chi) * L.u / (1 + L.chi);
        const Real dXy_uy = (LN.chi-LS.chi) * L.v / (1 + L.chi);
        const Real dXz_uz = (LB.chi-LF.chi) * L.w / (1 + L.chi);
        o(ix,iy,iz).p = L.p + fac*(dU+dV+dW - (dXx_ux+dXy_uy+dXz_uz));
      }
    }
    else
    {
      for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
      for(int iy=0; iy<FluidBlock::sizeY; ++iy)
      for(int ix=0; ix<FluidBlock::sizeX; ++ix)
      {
        const FluidElement &L =lab(ix,iy,iz);
        const FluidElement &LW=lab(ix-1,iy,iz), &LE=lab(ix+1,iy,iz);
        const FluidElement &LS=lab(ix,iy-1,iz), &LN=lab(ix,iy+1,iz);
        const FluidElement &LF=lab(ix,iy,iz-1), &LB=lab(ix,iy,iz+1);
        const Real dU = LE.u-LW.u, dV = LN.v-LS.v, dW = LB.w-LF.w;
        const Real dXx_ux = (LE.chi-LW.chi) * L.u / (1 + L.chi);
        const Real dXy_uy = (LN.chi-LS.chi) * L.v / (1 + L.chi);
        const Real dXz_uz = (LB.chi-LF.chi) * L.w / (1 + L.chi);
        const Real FADE = fade(info, ix, iy, iz);
        o(ix,iy,iz).p = FADE*( L.p + fac*(dU+dV+dW - (dXx_ux+dXy_uy+dXz_uz) ) );
      }
    }
  }
};

class OperatorGradP : public GenericLabOperator
{
 private:
  const double dt;
  const Real extent[3];

 public:
  OperatorGradP(double _dt,const Real ext[3]):dt(_dt),extent{ext[0],ext[1],ext[2]}
  {
    stencil = StencilInfo(-1,-1,-1, 2,2,2, false, 1, 4);
    stencil_start[0] = -1; stencil_start[1] = -1; stencil_start[2] = -1;
    stencil_end[0] = 2; stencil_end[1] = 2; stencil_end[2] = 2;
  }

  ~OperatorGradP() {}

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
      o(ix,iy,iz).u = (Uf+o(ix,iy,iz).chi*o(ix,iy,iz).tmpU)/(1+o(ix,iy,iz).chi);
      o(ix,iy,iz).v = (Vf+o(ix,iy,iz).chi*o(ix,iy,iz).tmpV)/(1+o(ix,iy,iz).chi);
      o(ix,iy,iz).w = (Wf+o(ix,iy,iz).chi*o(ix,iy,iz).tmpW)/(1+o(ix,iy,iz).chi);
    }
  }
};

template <typename Lab>
class CoordinatorPressure : public GenericCoordinator
{
 protected:
  PoissonSolver * pressureSolver;

 public:
  CoordinatorPressure(SimulationData & s) : GenericCoordinator(s)
  {
    if(sim.bUseFourierBC)
    pressureSolver = new PoissonSolverPeriodic<FluidGridMPI,StreamerDiv>(sim);
    else if (sim.bUseUnboundedBC)
    pressureSolver = new PoissonSolverUnbounded<FluidGridMPI,StreamerDiv>(sim);
    else
    pressureSolver = new PoissonSolverMixed<FluidGridMPI,StreamerDiv>(sim);
  }

  void operator()(const double dt)
  {
    check("pressure - start");
    const int nthreads = omp_get_max_threads();

    {
      std::vector<OperatorDivergence*> diff(nthreads, nullptr);
      #pragma omp parallel for schedule(static, 1)
      for(int i=0;i<nthreads;++i)
        diff[i] = new OperatorDivergence(dt, sim.fadeOutLength, sim.extent);

      compute<OperatorDivergence>(diff);
      for(int i=0; i<nthreads; i++) delete diff[i];
    }

    pressureSolver->solve();

    {
      //zero fields, going to contain Udef:
      #pragma omp parallel for schedule(static)
      for(unsigned i=0; i<vInfo.size(); i++) {
        const BlockInfo& info = vInfo[i];
        FluidBlock& b = *(FluidBlock*)info.ptrBlock;
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
    { //pressure correction dudt* = - grad P / rho
      std::vector<OperatorGradP*> diff(nthreads, nullptr);
      #pragma omp parallel for schedule(static, 1)
      for(int i=0;i<nthreads;++i) diff[i] = new OperatorGradP(dt, sim.extent);

      compute<OperatorGradP>(diff);
      for(int i=0; i<nthreads; i++) delete diff[i];
    }

    check("pressure - end");
  }

  std::string getName()
  {
    return "Pressure";
  }
};
#endif
