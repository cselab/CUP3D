//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#include "AdvectionDiffusion.h"

class KernelAdvectDiffuse
{
  private:
  const double dt;
  const double mu;
  const double lambda;
  const Real* const uInf;

  public:
  const std::array<int, 3> stencil_start = {-1, -1, -1};
  const std::array<int, 3> stencil_end = {2, 2, 2};
  const StencilInfo stencil = StencilInfo(-1,-1,-1, 2,2,2, false, 3, 1,2,3);

  KernelAdvectDiffuse(double _dt, double m, Real*const u, Real l) :
    dt(_dt), mu(m), lambda(l), uInf(u) { }

  template <typename Lab, typename BlockType>
  void operator()(Lab & lab, const BlockInfo& info, BlockType& o) const
  {
    const Real facA = -dt/(2*info.h_gridpoint);
    const Real facD = (mu/info.h_gridpoint) * (dt/info.h_gridpoint);
    for (int iz=0; iz<FluidBlock::sizeZ; ++iz)
    for (int iy=0; iy<FluidBlock::sizeY; ++iy)
    for (int ix=0; ix<FluidBlock::sizeX; ++ix)
    {
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
      o(ix,iy,iz).tmpU = L.u + facA*duA + facD*duD;
      o(ix,iy,iz).tmpV = L.v + facA*dvA + facD*dvD;
      o(ix,iy,iz).tmpW = L.w + facA*dwA + facD*dwD;
    }
  }
};

class KernelExplPenPart1
{
  private:
  const double dt;
  const double lambda;

  public:
  const std::array<int, 3> stencil_start = {-1, -1, -1};
  const std::array<int, 3> stencil_end = {2, 2, 2};
  const StencilInfo stencil = StencilInfo(-1,-1,-1, 2,2,2, false, 1, 4);

  KernelExplPenPart1(double _dt, Real l) : dt(_dt), lambda(l) { }

  template <typename Lab, typename BlockType>
  void operator()(Lab & lab, const BlockInfo& info, BlockType& o) const
  {
    const Real fac = - 0.5 * dt / info.h_gridpoint;
    for (int iz=0; iz<FluidBlock::sizeZ; ++iz)
    for (int iy=0; iy<FluidBlock::sizeY; ++iy)
    for (int ix=0; ix<FluidBlock::sizeX; ++ix) {
      const FluidElement &L =lab(ix,iy,iz);
      const FluidElement &LW=lab(ix-1,iy,iz), &LE=lab(ix+1,iy,iz);
      const FluidElement &LS=lab(ix,iy-1,iz), &LN=lab(ix,iy+1,iz);
      const FluidElement &LF=lab(ix,iy,iz-1), &LB=lab(ix,iy,iz+1);
      o(ix,iy,iz).u = L.tmpU - L.chi * L.u; // explicit penal part 1
      o(ix,iy,iz).v = L.tmpV - L.chi * L.v; // explicit penal part 1
      o(ix,iy,iz).w = L.tmpW - L.chi * L.w; // explicit penal part 1
      o(ix,iy,iz).tmpU = L.tmpU + fac*(LE.p-LW.p); // approximate ufuild
      o(ix,iy,iz).tmpV = L.tmpV + fac*(LN.p-LS.p); // approximate ufuild
      o(ix,iy,iz).tmpW = L.tmpW + fac*(LB.p-LF.p); // using old pressure
    }
  }
};

void AdvectionDiffusion::operator()(const double dt)
{
  sim.startProfiler("AdvDiff Kernel");
  {
    const KernelAdvectDiffuse K(dt, sim.nu, sim.uinf, sim.lambda);
    compute<KernelAdvectDiffuse>(K);
  }
  sim.stopProfiler();

  sim.startProfiler("AdvDiff copy");
  #if PENAL_TYPE==0
  {
    const KernelExplPenPart1 K(dt, sim.lambda);
    compute<KernelExplPenPart1>(K);
  }
  #else
    #pragma omp parallel for schedule(static)
    for(size_t i=0; i<vInfo.size(); i++)
    {
      FluidBlock& b = *(FluidBlock*) vInfo[i].ptrBlock;
      for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
      for(int iy=0; iy<FluidBlock::sizeY; ++iy)
      for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
        b(ix,iy,iz).u = b(ix,iy,iz).tmpU;
        b(ix,iy,iz).v = b(ix,iy,iz).tmpV;
        b(ix,iy,iz).w = b(ix,iy,iz).tmpW;
      }
    }
  #endif
  sim.stopProfiler();
  check("AdvectionDiffusion - end");
}
