//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch) and Christian Conti.
//

#ifndef CubismUP_3D_CoordinatorPressure_h
#define CubismUP_3D_CoordinatorPressure_h

#include "GenericOperator.h"
#include "GenericCoordinator.h"
#include "../obstacles/IF3D_ObstacleVector.h"
#ifdef _ACCFFT_
#ifdef CUP_UNBOUNDED_FFT
#include "../poisson/PoissonSolverScalarACC_freespace.h"
#else
#include "../poisson/PoissonSolverScalarACC.h"
#endif
typedef PoissonSolverScalarFFTW_ACC<FluidGridMPI, StreamerDiv> PressureSolver;
#else
#ifdef CUP_UNBOUNDED_FFT
#include "../poisson/PoissonSolverScalarFFTW_cyclicConvolution.h"
#else
#include "../poisson/PoissonSolverScalarFFTW.h"
#endif /* CUP_UNBOUNDED_FFT */
typedef PoissonSolverScalarFFTW_MPI<FluidGridMPI, StreamerDiv> PressureSolver;
#endif

class OperatorDivergence : public GenericLabOperator
{
 private:
  double dt;
  const Real BUF, iWidth = 1/BUF;
  const Real* const ext;
  PressureSolver*const solver;

  inline bool _is_touching(const BlockInfo& i) const
  {
    Real minP[2], maxP[2]; i.pos(minP, 0, 0, 0);
    i.pos(maxP, CUP_BLOCK_SIZE-1, CUP_BLOCK_SIZE-1, CUP_BLOCK_SIZE-1);
    const bool touchE = BUF>=ext[0]-maxP[0], touchW = BUF>=minP[0];
    const bool touchN = BUF>=ext[1]-maxP[1], touchS = BUF>=minP[1];
    #ifdef BC_PERIODICZ
    const bool touchF = false, touchB = false;
    #else
    const bool touchF = BUF>=ext[2]-maxP[2], touchB = BUF>=minP[2];
    #endif
    return touchN || touchE || touchS || touchW || touchF || touchB;
  }

 public:
  OperatorDivergence(double _dt, Real b, const Real* e, PressureSolver*const s)
   : dt(_dt), BUF(b), ext(e), solver(s)
  {
    stencil_start[0] = -1; stencil_start[1] = -1;  stencil_start[2] = -1;
    stencil_end[0] = 2;  stencil_end[1] = 2;  stencil_end[2] = 2;
    stencil = StencilInfo(-1,-1,-1, 2,2,2, false, 3, 1,2,3);
  }
  ~OperatorDivergence() {}

  template <typename Lab, typename BlockType>
  void operator()(Lab & lab, const BlockInfo& info, BlockType& o) const
  {
    const size_t offset = solver->_offset_ext(info);
    const Real h = info.h_gridpoint, fac = .5*h*h/dt;
    if( not _is_touching(info) )
    {
      for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
      for(int iy=0; iy<FluidBlock::sizeY; ++iy)
      for(int ix=0; ix<FluidBlock::sizeX; ++ix)
      {
        const Real dU = lab(ix+1,iy,iz).u - lab(ix-1,iy,iz).u;
        const Real dV = lab(ix,iy+1,iz).v - lab(ix,iy-1,iz).v;
        const Real dW = lab(ix,iy,iz+1).w - lab(ix,iy,iz-1).w;
        solver->_cub2fftw(offset, iz,iy,ix, lab(ix,iy,iz).p + fac*(dU+dV+dW));
      }
    }
    else
    {
      for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
      for(int iy=0; iy<FluidBlock::sizeY; ++iy)
      for(int ix=0; ix<FluidBlock::sizeX; ++ix)
      {
        Real p[3]; info.pos(p, ix, iy, iz);
        #ifdef BC_PERIODICZ
        const Real zt = 0, dzb = 0;
        #else
        const Real zt = std::max(Real(0), BUF -(ext[2]-p[2]) );
        const Real zb = std::max(Real(0), BUF - p[2] );
        #endif
        const Real yt = std::max(Real(0), BUF -(ext[1]-p[1]) );
        const Real yb = std::max(Real(0), BUF - p[1] );
        const Real xt = std::max(Real(0), BUF -(ext[0]-p[0]) );
        const Real xb = std::max(Real(0), BUF - p[0] );
        const Real dist = std::min(std::max({zt,zb,yt,yb,xt,xb}), BUF);
        const Real fade = 1 - std::pow(dist*iWidth, 2);
        const Real dU = lab(ix+1,iy,iz).u - lab(ix-1,iy,iz).u;
        const Real dV = lab(ix,iy+1,iz).v - lab(ix,iy-1,iz).v;
        const Real dW = lab(ix,iy,iz+1).w - lab(ix,iy,iz-1).w;
        const Real RHS = lab(ix,iy,iz).p + fac*(dU+dV+dW);
        solver->_cub2fftw(offset, iz,iy,ix, RHS*fade );
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
    for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
      // p contains the pressure correction after the Poisson solver
      o(ix,iy,iz).u += fac * (lab(ix+1,iy,iz).p - lab(ix-1,iy,iz).p);
      o(ix,iy,iz).v += fac * (lab(ix,iy+1,iz).p - lab(ix,iy-1,iz).p);
      o(ix,iy,iz).w += fac * (lab(ix,iy,iz+1).p - lab(ix,iy,iz-1).p);
    }
  }
};

template <typename Lab>
class CoordinatorPressure : public GenericCoordinator
{
 protected:
  PressureSolver pressureSolver;
  const Real * const buffer;
  IF3D_ObstacleVector** const obstacleVector;

 public:
  CoordinatorPressure(FluidGridMPI* g, IF3D_ObstacleVector** const ov,
    const Real * const b) :
    GenericCoordinator(g), pressureSolver(*g), buffer(b), obstacleVector(ov) {}

  void operator()(const double dt)
  {
    check("pressure - start");
    const int nthreads = omp_get_max_threads();

    #pragma omp parallel for schedule(static)
    for(unsigned i=0; i<vInfo.size(); i++) {
      const BlockInfo& info = vInfo[i];
      FluidBlock& b = *(FluidBlock*)info.ptrBlock;
      for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
      for(int iy=0; iy<FluidBlock::sizeY; ++iy)
      for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
        b(ix,iy,iz).tmpU = 0;
        b(ix,iy,iz).tmpV = 0;
        b(ix,iy,iz).tmpW = 0; //zero fields, going to contain Udef
      }
    }

    const Real h = grid->getH();
    const Real ext[3] = {
        h*grid->getBlocksPerDimension(0)*FluidBlock::sizeX,
        h*grid->getBlocksPerDimension(1)*FluidBlock::sizeY,
        h*grid->getBlocksPerDimension(2)*FluidBlock::sizeZ
    };


    {   //place onto p: ( div u^(t+1) - div u^* ) / dt
      //where i want div u^(t+1) to be equal to div udef
      std::vector<OperatorDivergence*> diff(nthreads, nullptr);
      #pragma omp parallel for schedule(static, 1)
      for(int i=0;i<nthreads;++i)
        diff[i] = new OperatorDivergence(dt, *buffer, ext, &pressureSolver);

      compute<OperatorDivergence>(diff);
      for(int i=0; i<nthreads; i++) delete diff[i];
    }

    pressureSolver.solve();

    { //pressure correction dudt* = - grad P / rho
      std::vector<OperatorGradP*> diff(nthreads, nullptr);
      #pragma omp parallel for schedule(static, 1)
      for(int i=0;i<nthreads;++i) diff[i] = new OperatorGradP(dt, ext);

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
