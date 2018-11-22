//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch) and Christian Conti.
//

#ifndef CubismUP_3D_CoordinatorAdvectDiffuse_h
#define CubismUP_3D_CoordinatorAdvectDiffuse_h

#include "GenericCoordinator.h"
#include "GenericOperator.h"
#include <cmath>

class OperatorAdvectDiffuse : public GenericLabOperator
{
  private:
  const double dt;
  const double mu;
  const Real* const uInf;

  public:
   OperatorAdvectDiffuse(const double _dt, double m, const Real* const u) :
   dt(_dt), mu(m), uInf(u) {
    stencil = StencilInfo(-1,-1,-1, 2,2,2, false, 3, 1,2,3);
    stencil_start[0] = -1; stencil_start[1] = -1; stencil_start[2] = -1;
    stencil_end[0] = 2; stencil_end[1] = 2;  stencil_end[2] = 2;
  }

  template <typename Lab, typename BlockType>
  void operator()(Lab & lab, const BlockInfo& info, BlockType& o)
  {
    const Real facA = -dt/(2*info.h_gridpoint);
    const Real facD = (mu/info.h_gridpoint) * (dt/info.h_gridpoint);

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

      const Real duD = LN.u+LS.u + LE.u+LW.u + LF.u+LB.u - L.u*6;
      const Real dvD = LN.v+LS.v + LE.v+LW.v + LF.v+LB.v - L.v*6;
      const Real dwD = LN.w+LS.w + LE.w+LW.w + LF.w+LB.w - L.w*6;

      const Real duA = u * dudx + v * dudy + w * dudz;
      const Real dvA = u * dvdx + v * dvdy + w * dvdz;
      const Real dwA = u * dwdx + v * dwdy + w * dwdz;

      o(ix,iy,iz).tmpU = L.u + facA * duA + facD * duD;
      o(ix,iy,iz).tmpV = L.v + facA * dvA + facD * dvD;
      o(ix,iy,iz).tmpW = L.w + facA * dwA + facD * dwD;
    }
  }
};


template <typename Lab>
class CoordinatorAdvectDiffuse : public GenericCoordinator
{
protected:
  const Real MU;
  const Real* const uInf;

public:
  CoordinatorAdvectDiffuse(const Real _mu, const Real* const _uInf,
    FluidGridMPI * _grid) : GenericCoordinator(_grid), MU(_mu), uInf(_uInf) { }

  ~CoordinatorAdvectDiffuse() { }

  void operator()(const double dt)
  {
    check("AdvectDiffuse - start");
    const int nthreads = omp_get_max_threads();

    {
      using advection = OperatorAdvectDiffuse;
      std::vector<advection*> adv1(nthreads, nullptr);
      #pragma omp parallel for schedule(static, 1)
      for(int i=0;i<nthreads;++i) adv1[i] = new advection(dt, MU, uInf);

      compute(adv1);

      for(int i=0; i<nthreads; i++) delete adv1[i];
    }

    #pragma omp parallel for schedule(static)
    for(size_t i=0; i<vInfo.size(); i++) {
      const BlockInfo& info = vInfo[i];
      FluidBlock& b = *(FluidBlock*)info.ptrBlock;
      for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
      for(int iy=0; iy<FluidBlock::sizeY; ++iy)
      for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
        b(ix,iy,iz).u = b(ix,iy,iz).tmpU;
        b(ix,iy,iz).v = b(ix,iy,iz).tmpV;
        b(ix,iy,iz).w = b(ix,iy,iz).tmpW;
      }
    }

    check("AdvectDiffuse - end");
  }

  std::string getName()
  {
    return "AdvectDiffuse";
  }
};
#endif
