//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch) and Christian Conti.
//

#ifndef CubismUP_3D_CoordinatorAdvection_h
#define CubismUP_3D_CoordinatorAdvection_h

#include "GenericCoordinator.h"
#include "GenericOperator.h"
#include <cmath>

class OperatorAdvection : public GenericLabOperator
{
  private:
  const double dt;
  const Real* const uInf;

  public:
   OperatorAdvection(const double _dt, const Real* const u) : dt(_dt), uInf(u)
  {
    stencil = StencilInfo(-1,-1,-1, 2,2,2, false, 3, 1,2,3);
    stencil_start[0] = -1; stencil_start[1] = -1; stencil_start[2] = -1;
    stencil_end[0] = 2; stencil_end[1] = 2;  stencil_end[2] = 2;
  }

  template <typename Lab, typename BlockType>
  void operator()(Lab & lab, const BlockInfo& info, BlockType& o)
  {
    #ifndef _RK2_
    const Real fac = -dt/(2.*info.h_gridpoint);
    #else //perform half step
    const Real fac = -dt/(4.*info.h_gridpoint);
    #endif

    for (int iz=0; iz<FluidBlock::sizeZ; ++iz)
    for (int iy=0; iy<FluidBlock::sizeY; ++iy)
    for (int ix=0; ix<FluidBlock::sizeX; ++ix) {
      const FluidElement &L =lab(ix,iy,iz);
      const FluidElement &LW=lab(ix-1,iy,iz), &LE=lab(ix+1,iy,iz);
      const FluidElement &LS=lab(ix,iy-1,iz), &LN=lab(ix,iy+1,iz);
      const FluidElement &LF=lab(ix,iy,iz-1), &LB=lab(ix,iy,iz+1);

      const Real dudx= LE.u-LW.u, dvdx= LE.v-LW.v, dwdx= LE.w-LW.w;
      const Real dudy= LN.u-LS.u, dvdy= LN.v-LS.v, dwdy= LN.w-LS.w;
      const Real dudz= LB.u-LF.u, dvdz= LB.v-LF.v, dwdz= LB.w-LF.w;
      const Real u = L.u+uInf[0], v = L.v+uInf[1], w = L.w+uInf[2];
      o(ix,iy,iz).tmpU = L.u +fac*(u * dudx + v * dudy + w * dudz);
      o(ix,iy,iz).tmpV = L.v +fac*(u * dvdx + v * dvdy + w * dvdz);
      o(ix,iy,iz).tmpW = L.w +fac*(u * dwdx + v * dwdy + w * dwdz);
    }
  }
};

#ifdef _RK2_
class OperatorAdvectionStage2 : public GenericLabOperator
{
  private:
  const double dt;
  const Real* const uInf;

  public:
  OperatorAdvectionStage2(const double _dt, const Real* const _uInf)
  : dt(_dt), uInf(_uInf)
  {
    stencil = StencilInfo(-1,-1,-1, 2,2,2, false, 3, 5,6,7);
    stencil_start[0] = -1; stencil_start[1] = -1; stencil_start[2] = -1;
    stencil_end[0] = 2;  stencil_end[1] = 2;  stencil_end[2] = 2;
  }

  ~OperatorAdvectionStage2() {}

  template <typename Lab, typename BlockType>
  void operator()(Lab & lab, const BlockInfo& info, BlockType& o) const
  {
    const Real fac = -dt/(2.*info.h_gridpoint);

    for (int iz=0; iz<FluidBlock::sizeZ; ++iz)
    for (int iy=0; iy<FluidBlock::sizeY; ++iy)
    for (int ix=0; ix<FluidBlock::sizeX; ++ix) {
      const FluidElement &L =lab(ix,iy,iz);
      const FluidElement &LW=lab(ix-1,iy,iz), &LE=lab(ix+1,iy,iz);
      const FluidElement &LS=lab(ix,iy-1,iz), &LN=lab(ix,iy+1,iz);
      const FluidElement &LF=lab(ix,iy,iz-1), &LB=lab(ix,iy,iz+1);

      const Real dudx=LE.tmpU-LW.tmpU, dvdx=LE.tmpV-LW.tmpV, dwdx=LE.tmpW-LW.tmpW;
      const Real dudy=LN.tmpU-LS.tmpU, dvdy=LN.tmpV-LS.tmpV, dwdy=LN.tmpW-LS.tmpW;
      const Real dudz=LB.tmpU-LF.tmpU, dvdz=LB.tmpV-LF.tmpV, dwdz=LB.tmpW-LF.tmpW;
      const Real u = L.u+uInf[0], v = L.v+uInf[1], w = L.w+uInf[2];

      o(ix,iy,iz).u = L.u +fac*(u * dudx + v * dudy + w * dudz);
      o(ix,iy,iz).v = L.v +fac*(u * dvdx + v * dvdy + w * dvdz);
      o(ix,iy,iz).w = L.w +fac*(u * dwdx + v * dwdy + w * dwdz);
    }
  }
};
#endif // _RK2_

class OperatorAdvectionUpwind3rdOrder : public GenericLabOperator
{
  private:
  const double dt;
  const Real* const uInf;

  public:
  OperatorAdvectionUpwind3rdOrder(const double _dt, const Real* const _uInf)
  : dt(_dt), uInf(_uInf)
  {
    stencil = StencilInfo(-2,-2,-2, 3,3,3, false, 3, 1,2,3);
    stencil_start[0] = -2; stencil_start[1] = -2; stencil_start[2] = -2;
    stencil_end[0] = 3;  stencil_end[1] = 3;  stencil_end[2] = 3;
  }

  ~OperatorAdvectionUpwind3rdOrder() {}

  template <typename Lab, typename BlockType>
  void operator()(Lab & lab, const BlockInfo& info, BlockType& o) const
  {
    #ifndef _RK2_
    const Real factor = -dt/(6.*info.h_gridpoint);
    #else //perform half step
    const Real factor = -dt/(12.*info.h_gridpoint);
    #endif

    for (int iz=0; iz<FluidBlock::sizeZ; ++iz)
    for (int iy=0; iy<FluidBlock::sizeY; ++iy)
    for (int ix=0; ix<FluidBlock::sizeX; ++ix)
    {
      const FluidElement& L   = lab(ix  ,iy  ,iz  );
      const FluidElement &LW =lab(ix-1,iy,iz), &LE =lab(ix+1,iy,iz);
      const FluidElement &LS =lab(ix,iy-1,iz), &LN =lab(ix,iy+1,iz);
      const FluidElement &LF =lab(ix,iy,iz-1), &LB =lab(ix,iy,iz+1);
      const FluidElement &LW2=lab(ix-2,iy,iz), &LE2=lab(ix+2,iy,iz);
      const FluidElement &LS2=lab(ix,iy-2,iz), &LN2=lab(ix,iy+2,iz);
      const FluidElement &LF2=lab(ix,iy,iz-2), &LB2=lab(ix,iy,iz+2);

      const Real u3 = 3*L.u, v3 = 3*L.v, w3 = 3*L.w;
      const Real u = L.u+uInf[0], v = L.v+uInf[1], w = L.w+uInf[2];

      const Real dudx = u>0 ?         2*LE.u +u3 -6*LW.u +LW2.u
                            : -LE2.u +6*LE.u -u3 -2*LW.u;
      const Real dvdx = u>0 ?         2*LE.v +v3 -6*LW.v +LW2.v
                            : -LE2.v +6*LE.v -v3 -2*LW.v;
      const Real dwdx = u>0 ?         2*LE.w +w3 -6*LW.w +LW2.w
                            : -LE2.w +6*LE.w -w3 -2*LW.w;

      const Real dudy = v>0 ?         2*LN.u +u3 -6*LS.u +LS2.u
                            : -LN2.u +6*LN.u -u3 -2*LS.u;
      const Real dvdy = v>0 ?         2*LN.v +v3 -6*LS.v +LS2.v
                            : -LN2.v +6*LN.v -v3 -2*LS.v;
      const Real dwdy = v>0 ?         2*LN.w +w3 -6*LS.w +LS2.w
                            : -LN2.w +6*LN.w -w3 -2*LS.w;

      const Real dudz = w>0 ?         2*LB.u +u3 -6*LF.u +LF2.u
                            : -LB2.u +6*LB.u -u3 -2*LF.u;
      const Real dvdz = w>0 ?         2*LB.v +v3 -6*LF.v +LF2.v
                            : -LB2.v +6*LB.v -v3 -2*LF.v;
      const Real dwdz = w>0 ?         2*LB.w +w3 -6*LF.w +LF2.w
                            : -LB2.w +6*LB.w -w3 -2*LF.w;

      o(ix,iy,iz).tmpU = L.u + factor*(u * dudx + v * dudy + w * dudz);
      o(ix,iy,iz).tmpV = L.v + factor*(u * dvdx + v * dvdy + w * dvdz);
      o(ix,iy,iz).tmpW = L.w + factor*(u * dwdx + v * dwdy + w * dwdz);
    }
  }
};

#ifdef _RK2_
class OperatorAdvectionUpwind3rdOrderStage2 : public GenericLabOperator
{
  private:
  const double dt;
  const Real* const uInf;

  public:
  OperatorAdvectionUpwind3rdOrderStage2(const double dt, const Real* const uInf)
  : dt(_dt), uInf(_uInf)
  {
    stencil = StencilInfo(-2,-2,-2, 3,3,3, false, 3, 5,6,7);
    stencil_start[0] = -2; stencil_start[1] = -2; stencil_start[2] = -2;
    stencil_end[0] = 3; stencil_end[1] = 3; stencil_end[2] = 3;
  }

  ~OperatorAdvectionUpwind3rdOrderStage2() {}

  template <typename Lab, typename BlockType>
  void operator()(Lab & lab, const BlockInfo& info, BlockType& o) const
  {
    const Real factor = -dt/(6.*info.h_gridpoint);

    for (int iz=0; iz<FluidBlock::sizeZ; ++iz)
    for (int iy=0; iy<FluidBlock::sizeY; ++iy)
    for (int ix=0; ix<FluidBlock::sizeX; ++ix)
    {
      const FluidElement& L  = lab(ix,iy,iz);
      const FluidElement &LW =lab(ix-1,iy,iz), &LE =lab(ix+1,iy,iz);
      const FluidElement &LS =lab(ix,iy-1,iz), &LN =lab(ix,iy+1,iz);
      const FluidElement &LF =lab(ix,iy,iz-1), &LB =lab(ix,iy,iz+1);
      const FluidElement &LW2=lab(ix-2,iy,iz), &LE2=lab(ix+2,iy,iz);
      const FluidElement &LS2=lab(ix,iy-2,iz), &LN2=lab(ix,iy+2,iz);
      const FluidElement &LF2=lab(ix,iy,iz-2), &LB2=lab(ix,iy,iz+2);
      const Real u3 = 3*L.tmpU, v3 = 3*L.tmpV, w3 = 3*L.tmpW;
      const Real u = L.u+uInf[0], v = L.v+uInf[1], w = L.w+uInf[2];

      const Real dudx = u>0 ?            2*LE.tmpU +u3 -6*LW.tmpU +LW2.tmpU
                            : -LE2.tmpU +6*LE.tmpU -u3 -2*LW.tmpU;
      const Real dvdx = v>0 ?            2*LE.tmpV +v3 -6*LW.tmpV +LW2.tmpV
                            : -LE2.tmpV +6*LE.tmpV -v3 -2*LW.tmpV;
      const Real dwdx = w>0 ?            2*LE.tmpW +w3 -6*LW.tmpW +LW2.tmpW
                            : -LE2.tmpW +6*LE.tmpW -w3 -2*LW.tmpW;

      const Real dudy = u>0 ?            2*LN.tmpU +u3 -6*LS.tmpU +LS2.tmpU
                            : -LN2.tmpU +6*LN.tmpU -u3 -2*LS.tmpU;
      const Real dvdy = v>0 ?            2*LN.tmpV +v3 -6*LS.tmpV +LS2.tmpV
                            : -LN2.tmpV +6*LN.tmpV -v3 -2*LS.tmpV;
      const Real dwdy = w>0 ?            2*LN.tmpW +w3 -6*LS.tmpW +LS2.tmpW
                            : -LN2.tmpW +6*LN.tmpW -w3 -2*LS.tmpW;

      const Real dudz = u>0 ?            2*LB.tmpU +u3 -6*LF.tmpU +LF2.tmpU
                            : -LB2.tmpU +6*LB.tmpU -u3 -2*LF.tmpU;
      const Real dvdz = v>0 ?            2*LB.tmpV +v3 -6*LF.tmpV +LF2.tmpV
                            : -LB2.tmpV +6*LB.tmpV -v3 -2*LF.tmpV;
      const Real dwdz = w>0 ?            2*LB.tmpW +w3 -6*LF.tmpW +LF2.tmpW
                            : -LB2.tmpW +6*LB.tmpW -w3 -2*LF.tmpW;

      o(ix,iy,iz).u = L.u + factor*(u * dudx + v * dudy + w * dudz);
      o(ix,iy,iz).v = L.v + factor*(u * dvdx + v * dvdy + w * dvdz);
      o(ix,iy,iz).w = L.w + factor*(u * dwdx + v * dwdy + w * dwdz);
    }
  }
};
#endif // _RK2_

template <typename Lab>
class CoordinatorAdvection : public GenericCoordinator
{
protected:
  const Real* const uInf;

public:
  CoordinatorAdvection(const Real* const _uInf, FluidGridMPI * _grid)
  : GenericCoordinator(_grid), uInf(_uInf)
  { }

  ~CoordinatorAdvection()
  { }

  void operator()(const double dt)
  {
    check("advection - start");
    const int nthreads = omp_get_max_threads();

    {
      //using advection = OperatorAdvectionUpwind3rdOrder;
      using advection = OperatorAdvection;
      std::vector<advection*> adv1(nthreads, nullptr);
      for(int i=0;i<nthreads;++i) adv1[i] = new advection(dt, uInf);

      compute(adv1);

      for(int i=0; i<nthreads; i++) delete adv1[i];
    }
#ifdef _RK2_
    {
      //using advection = OperatorAdvectionUpwind3rdOrderStage2;
      using advection = OperatorAdvectionStage2;
      std::vector<advection*> adv2(nthreads, nullptr);
      for(int i=0;i<nthreads;++i) adv2[i] = new advection(dt, uInf);

      compute(adv2);

      for(int i=0; i<nthreads; i++) delete adv2[i];
    }
#endif

    check("advection - end");
  }

  std::string getName()
  {
    return "Advection";
  }
};
#endif
