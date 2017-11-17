//
//  CubismUP_3D
//
//  Written by Guido Novati ( novatig@ethz.ch ).
//  This file started as an extension of code written by Christian Conti
//  Copyright (c) 2017 ETHZ. All rights reserved.
//

#ifndef CubismUP_3D_CoordinatorDiffusion_h
#define CubismUP_3D_CoordinatorDiffusion_h

#include "GenericCoordinator.h"
#include "GenericOperator.h"

class OperatorDiffusion_toField : public GenericLabOperator
{
private:
  const double mu;
  const double dt;

public:
  OperatorDiffusion_toField(double _dt, double m) : mu(m), dt(_dt)
  {
    stencil = StencilInfo(-1,-1,-1, 2,2,2, false, 3, 5,6,7);
    stencil_start[0] = stencil_start[1] = stencil_start[2] = -1;
    stencil_end[0] = stencil_end[1] = stencil_end[2] = 2;
  }

  ~OperatorDiffusion_toField() {}

  template <typename Lab, typename BlockType>
  void operator()(Lab & lab, const BlockInfo& info, BlockType& o) const
  {
#ifdef _RK2_
    const Real fac = 0.5 * (mu/info.h_gridpoint) * (dt/info.h_gridpoint);
#else
    const Real fac = (mu/info.h_gridpoint) * (dt/info.h_gridpoint);
#endif
    for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
    for(int iy=0; iy<FluidBlock::sizeY; ++iy)
    for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
      const FluidElement &L =lab(ix,iy,iz);
      const FluidElement &LW=lab(ix-1,iy,iz), &LE=lab(ix+1,iy,iz);
      const FluidElement &LS=lab(ix,iy-1,iz), &LN=lab(ix,iy+1,iz);
      const FluidElement &LF=lab(ix,iy,iz-1), &LB=lab(ix,iy,iz+1);
      #ifdef _RK2_
      const Real u = L.u, v = L.v, w = L.w;
      #else
      const Real u = L.tmpU, v = L.tmpV, w = L.tmpW;
      #endif
      o(ix,iy,iz).u = u + fac*(LN.tmpU+LS.tmpU + LE.tmpU+LW.tmpU
                             + LF.tmpU+LB.tmpU - L.tmpU*6);
      o(ix,iy,iz).v = v + fac*(LN.tmpV+LS.tmpV + LE.tmpV+LW.tmpV
                             + LF.tmpV+LB.tmpV - L.tmpV*6);
      o(ix,iy,iz).w = w + fac*(LN.tmpW+LS.tmpW + LE.tmpW+LW.tmpW
                             + LF.tmpW+LB.tmpW - L.tmpW*6);
    }
  }
};

#ifdef _RK2_
class OperatorDiffusion_toTemp : public GenericLabOperator
{
private:
  const double mu;
  const double dt;

public:
  OperatorDiffusion_toTemp(double _dt, double m) : mu(m), dt(_dt)
  {
    stencil = StencilInfo(-1,-1,-1, 2,2,2, false, 3, 1,2,3);
    stencil_start[0] = stencil_start[1] = stencil_start[2] = -1;
    stencil_end[0] = stencil_end[1] = stencil_end[2] = 2;
  }

  ~OperatorDiffusion_toTemp() {}

  template <typename Lab, typename BlockType>
  void operator()(Lab & lab, const BlockInfo& info, BlockType& o) const
  {
    const Real fac = 0.5 * (mu/info.h_gridpoint) * (dt/info.h_gridpoint);
    for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
    for(int iy=0; iy<FluidBlock::sizeY; ++iy)
    for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
      const FluidElement &L =lab(ix,iy,iz);
      const FluidElement &LW=lab(ix-1,iy,iz), &LE=lab(ix+1,iy,iz);
      const FluidElement &LS=lab(ix,iy-1,iz), &LN=lab(ix,iy+1,iz);
      const FluidElement &LF=lab(ix,iy,iz-1), &LB=lab(ix,iy,iz+1);

      o(ix,iy,iz).tmpU = L.u +fac*(LN.u+LS.u + LE.u+LW.u + LF.u+LB.u - L.u*6);
      o(ix,iy,iz).tmpV = L.v +fac*(LN.v+LS.v + LE.v+LW.v + LF.v+LB.v - L.v*6);
      o(ix,iy,iz).tmpW = L.w +fac*(LN.w+LS.w + LE.w+LW.w + LF.w+LB.w - L.w*6);
    }
  }
};
#endif // _RK2_

template <typename Lab>
class CoordinatorDiffusion : public GenericCoordinator
{
protected:
  const Real coeff;

public:
  CoordinatorDiffusion(const Real coeff, FluidGridMPI * grid) : GenericCoordinator(grid), coeff(coeff)
  {
  }

  void operator()(const Real dt)
  {
    check("diffusion - start");
    const int nthreads = omp_get_max_threads();

#ifdef _RK2_
    {
      vector<OperatorDiffusion_toTemp*> diff1(nthreads, nullptr);
      for(int i=0;i<nthreads;++i)
        diff1[i] = new OperatorDiffusion_toTemp(dt, coeff);

      compute(diff1);

      for(int i=0; i<nthreads; i++) delete diff1[i];
    }
#endif
    {
      vector<OperatorDiffusion_toField*> diff2(nthreads, nullptr);
      for(int i=0;i<nthreads;++i)
        diff2[i] = new OperatorDiffusion_toField(dt, coeff);

      compute(diff2);

      for(int i=0; i<nthreads; i++) delete diff2[i];
    }

    check("diffusion - end");
  }

  string getName()
  {
    return "Diffusion";
  }
};

#endif
