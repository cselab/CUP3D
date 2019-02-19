//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch) and Christian Conti.
//

#ifndef CubismUP_3D_CoordinatorPressureRHS_h
#define CubismUP_3D_CoordinatorPressureRHS_h

#include "PenalizationObstacleVisitor.h"
#include "../obstacles/IF3D_ObstacleVector.h"

//#include "../poisson/PoissonSolverScalarACC_freespace.h"
//#include "../poisson/PoissonSolverScalarACC.h"

#include "../poisson/PoissonSolverUnbounded.h"
#include "../poisson/PoissonSolverPeriodic.h"
#include "../poisson/PoissonSolverMixed.h"

class OperatorPressureRHS : public GenericLabOperator
{
 private:
  double dt;
  const Real ext[3], fadeLen[3], iFade[3];
  static constexpr Real EPS = std::numeric_limits<Real>::epsilon();
  PoissonSolver * const solver;

  inline bool _is_touching(const BlockInfo& i) const
  {
    Real maxP[3], minP[3]; i.pos(minP, 0, 0, 0);
    i.pos(maxP, CUP_BLOCK_SIZE-1, CUP_BLOCK_SIZE-1, CUP_BLOCK_SIZE-1);
    const bool touchW= fadeLen[0]>=minP[0], touchE= fadeLen[0]>=ext[0]-maxP[0];
    const bool touchS= fadeLen[1]>=minP[1], touchN= fadeLen[1]>=ext[1]-maxP[1];
    const bool touchB= fadeLen[2]>=minP[2], touchF= fadeLen[2]>=ext[2]-maxP[2];
    return touchN || touchE || touchS || touchW || touchF || touchB;
  }

  inline Real fade(const BlockInfo&i, const int x,const int y,const int z) const
  {
    Real p[3]; i.pos(p, x, y, z);
    const Real zt = iFade[2] * std::max(Real(0), fadeLen[2] -(ext[2]-p[2]) );
    const Real zb = iFade[2] * std::max(Real(0), fadeLen[2] - p[2] );
    const Real yt = iFade[1] * std::max(Real(0), fadeLen[1] -(ext[1]-p[1]) );
    const Real yb = iFade[1] * std::max(Real(0), fadeLen[1] - p[1] );
    const Real xt = iFade[0] * std::max(Real(0), fadeLen[0] -(ext[0]-p[0]) );
    const Real xb = iFade[0] * std::max(Real(0), fadeLen[0] - p[0] );
    return 1-std::pow(std::min( std::max({zt,zb,yt,yb,xt,xb}), (Real)1), 2);
  }

  inline Real RHS(Lab&l, const int x,const int y,const int z,const Real F) const
  {
    const FluidElement & L  = l(x,  y,  z);
    const FluidElement & LW = l(x-1,y,  z  ), & LE = l(x+1,y,  z  );
    const FluidElement & LS = l(x,  y-1,z  ), & LN = l(x,  y+1,z  );
    const FluidElement & LF = l(x,  y,  z-1), & LB = l(x,  y,  z+1);
    const Real divUs = LE.tmpU-LW.tmpU + LN.tmpV-LS.tmpV + LB.tmpW-LF.tmpW;
    const Real divUt = LE.u-LW.u + LN.v-LS.v + LB.w-LF.w, X = L.chi;
    #ifndef ZEROGRADCHI
      const Real dXx_dux = (LE.chi-LW.chi)*( L.tmpU-L.u + F * (LE.p-LW.p) );
      const Real dXy_duy = (LN.chi-LS.chi)*( L.tmpV-L.v + F * (LN.p-LS.p) );
      const Real dXz_duz = (LB.chi-LF.chi)*( L.tmpW-L.w + F * (LB.p-LF.p) );
      return divUt -X*X*divUs + (dXx_dux+dXy_duy+dXz_duz) / (1 + X);
    #else
      return divUt -X*X*divUs;
    #endif
  }

 public:
  OperatorPressureRHS(double _dt, const Real buf[3], const Real extent[3],
   PoissonSolver * ps) : dt(_dt), ext{extent[0],extent[1],extent[2]},
   fadeLen{buf[0],buf[1],buf[2]}, iFade{1/(buf[0]+EPS), 1/(buf[1]+EPS),
     1/(buf[2]+EPS)}, solver(ps)
  {
    stencil_start[0] = -1; stencil_start[1] = -1;  stencil_start[2] = -1;
    stencil_end[0] = 2;  stencil_end[1] = 2;  stencil_end[2] = 2;
    stencil = StencilInfo(-1,-1,-1, 2,2,2, false, 8, 0,1,2,3,4,5,6,7);
  }
  ~OperatorPressureRHS() {}

  template <typename Lab, typename BlockType>
  void operator()(Lab & lab, const BlockInfo& info, BlockType& o) const
  {
    const Real h = info.h_gridpoint, fac = .5*h*h/dt, pFac = .5*dt/h;
    const size_t offset = solver->_offset_ext(info);
    if( not _is_touching(info) )
    {
      for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
      for(int iy=0; iy<FluidBlock::sizeY; ++iy)
      for(int ix=0; ix<FluidBlock::sizeX; ++ix)
        solver->_cub2fftw(offset, iz,iy,ix, fac * RHS(lab, ix,iy,iz, pFac) );
    }
    else
    {
      for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
      for(int iy=0; iy<FluidBlock::sizeY; ++iy)
      for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
        const Real RHS_ = fade(info, ix,iy,iz) * RHS(lab, ix,iy,iz, pFac);
        solver->_cub2fftw(offset, iz,iy,ix, fac * RHS_);
      }
    }
  }
};

template <typename Lab>
class CoordinatorPressureRHS : public GenericCoordinator
{
 public:
  CoordinatorPressureRHS(SimulationData & s) : GenericCoordinator(s) {}

  void operator()(const double dt)
  {
    const int nthreads = omp_get_max_threads();
    {
      //zero fields, going to contain Udef:
      #pragma omp parallel for schedule(static)
      for(unsigned i=0; i<vInfo.size(); i++)
      {
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

    //place onto p: ( div u^(t+1) - div u^* ) / dt
    //where i want div u^(t+1) to be equal to div udef
    std::vector<OperatorPressureRHS*> diff(nthreads, nullptr);
    #pragma omp parallel for schedule(static, 1)
    for(int i=0; i<nthreads; ++i)
      diff[i] = new OperatorPressureRHS(dt, sim.fadeOutLengthPRHS, sim.extent, sim.pressureSolver);

    compute<OperatorPressureRHS>(diff);
    for(int i=0; i<nthreads; i++) delete diff[i];

    check("pressure rhs - end");
  }

  std::string getName() { return "PressureRHS"; }
};
#endif
