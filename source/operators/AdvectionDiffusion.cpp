//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#include "operators/AdvectionDiffusion.h"

CubismUP_3D_NAMESPACE_BEGIN
using namespace cubism;

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

class KernelAdvectDiffuse_HighOrder
{
  #define FOURTH_ORD_ADV
  private:
  const double dt;
  const double mu;
  const double lambda;
  const Real* const uInf;

  public:
  const std::array<int, 3> stencil_start = {-2,-2,-2}, stencil_end = {3, 3, 3};
  const StencilInfo stencil = StencilInfo(-2,-2,-2, 3,3,3, false, 3, 1,2,3);

  KernelAdvectDiffuse_HighOrder(double _dt, double m, Real*const u, Real l) :
    dt(_dt), mu(m), lambda(l), uInf(u) { }

  template <typename Lab, typename BlockType>
  void operator()(Lab & lab, const BlockInfo& info, BlockType& o) const
  {
    #ifdef FOURTH_ORD_ADV
      const Real facA = -dt/(12*info.h_gridpoint);
    #else
      const Real facA = -dt/(6*info.h_gridpoint);
    #endif
    const Real facD = (mu/info.h_gridpoint) * (dt/info.h_gridpoint) / 12;
    for (int iz=0; iz<FluidBlock::sizeZ; ++iz)
    for (int iy=0; iy<FluidBlock::sizeY; ++iy)
    for (int ix=0; ix<FluidBlock::sizeX; ++ix)
    {
      const FluidElement &L =lab(ix,iy,iz);
      const FluidElement &LW =lab(ix-1,iy,iz), &LE =lab(ix+1,iy,iz);
      const FluidElement &LS =lab(ix,iy-1,iz), &LN =lab(ix,iy+1,iz);
      const FluidElement &LF =lab(ix,iy,iz-1), &LB =lab(ix,iy,iz+1);
      const FluidElement &LW2=lab(ix-2,iy,iz), &LE2=lab(ix+2,iy,iz);
      const FluidElement &LS2=lab(ix,iy-2,iz), &LN2=lab(ix,iy+2,iz);
      const FluidElement &LF2=lab(ix,iy,iz-2), &LB2=lab(ix,iy,iz+2);
      const Real u = L.u+uInf[0], v = L.v+uInf[1], w = L.w+uInf[2];
      #ifdef FOURTH_ORD_ADV
        const Real dudx =  -LE2.u + 8*LE.u - 8*LW.u + LW2.u;
        const Real dvdx =  -LE2.v + 8*LE.v - 8*LW.v + LW2.v;
        const Real dwdx =  -LE2.w + 8*LE.w - 8*LW.w + LW2.w;
        const Real dudy =  -LN2.u + 8*LN.u - 8*LS.u + LS2.u;
        const Real dvdy =  -LN2.v + 8*LN.v - 8*LS.v + LS2.v;
        const Real dwdy =  -LN2.w + 8*LN.w - 8*LS.w + LS2.w;
        const Real dudz =  -LB2.u + 8*LB.u - 8*LF.u + LF2.u;
        const Real dvdz =  -LB2.v + 8*LB.v - 8*LF.v + LF2.v;
        const Real dwdz =  -LB2.w + 8*LB.w - 8*LF.w + LF2.w;
      #else // third order upwind
        const Real dudx = u>0 ?         2*LE.u +3*L.u -6*LW.u +LW2.u
                              : -LE2.u +6*LE.u -3*L.u -2*LW.u;
        const Real dvdx = u>0 ?         2*LE.v +3*L.v -6*LW.v +LW2.v
                              : -LE2.v +6*LE.v -3*L.v -2*LW.v;
        const Real dwdx = u>0 ?         2*LE.w +3*L.w -6*LW.w +LW2.w
                              : -LE2.w +6*LE.w -3*L.w -2*LW.w;
        const Real dudy = v>0 ?         2*LN.u +3*L.u -6*LS.u +LS2.u
                              : -LN2.u +6*LN.u -3*L.u -2*LS.u;
        const Real dvdy = v>0 ?         2*LN.v +3*L.v -6*LS.v +LS2.v
                              : -LN2.v +6*LN.v -3*L.v -2*LS.v;
        const Real dwdy = v>0 ?         2*LN.w +3*L.w -6*LS.w +LS2.w
                              : -LN2.w +6*LN.w -3*L.w -2*LS.w;
        const Real dudz = w>0 ?         2*LB.u +3*L.u -6*LF.u +LF2.u
                              : -LB2.u +6*LB.u -3*L.u -2*LF.u;
        const Real dvdz = w>0 ?         2*LB.v +3*L.v -6*LF.v +LF2.v
                              : -LB2.v +6*LB.v -3*L.v -2*LF.v;
        const Real dwdz = w>0 ?         2*LB.w +3*L.w -6*LF.w +LF2.w
                              : -LB2.w +6*LB.w -3*L.w -2*LF.w;
      #endif
      const Real dudx2 = 16*(LW.u+LE.u) -(LW2.u+LE2.u) - 30*L.u;
      const Real dvdx2 = 16*(LW.v+LE.v) -(LW2.v+LE2.v) - 30*L.v;
      const Real dwdx2 = 16*(LW.w+LE.w) -(LW2.w+LE2.w) - 30*L.w;
      const Real dudy2 = 16*(LS.u+LN.u) -(LS2.u+LN2.u) - 30*L.u;
      const Real dvdy2 = 16*(LS.v+LN.v) -(LS2.v+LN2.v) - 30*L.v;
      const Real dwdy2 = 16*(LS.w+LN.w) -(LS2.w+LN2.w) - 30*L.w;
      const Real dudz2 = 16*(LF.u+LB.u) -(LF2.u+LB2.u) - 30*L.u;
      const Real dvdz2 = 16*(LF.v+LB.v) -(LF2.v+LB2.v) - 30*L.v;
      const Real dwdz2 = 16*(LF.w+LB.w) -(LF2.w+LB2.w) - 30*L.w;
      const Real duA = u * dudx + v * dudy + w * dudz;
      const Real dvA = u * dvdx + v * dvdy + w * dvdz;
      const Real dwA = u * dwdx + v * dwdy + w * dwdz;
      o(ix,iy,iz).tmpU = L.u + facA*duA + facD*(dudx2 + dudy2 + dudz2);
      o(ix,iy,iz).tmpV = L.v + facA*dvA + facD*(dvdx2 + dvdy2 + dvdz2);
      o(ix,iy,iz).tmpW = L.w + facA*dwA + facD*(dwdx2 + dwdy2 + dwdz2);
    }
  }
};

class KernelAdvectDiffuse_nonUniform
{
  private:
  const Real dt, mu;
  const double lambda;
  const Real* const uInf;

  public:
  const std::array<int, 3> stencil_start = {-1, -1, -1};
  const std::array<int, 3> stencil_end = {2, 2, 2};
  const StencilInfo stencil = StencilInfo(-1,-1,-1, 2,2,2, false, 3, 1,2,3);

  KernelAdvectDiffuse_nonUniform(double _dt, double m, Real*const u, Real l) :
    dt(_dt), mu(m), lambda(l), uInf(u) { }

  template <typename Lab, typename BlockType>
  void operator()(Lab & lab, const BlockInfo& info, BlockType& o) const
  {
    // FD coefficients for first and second derivative
    const BlkCoeffX & c1x = o.fd_cx.first, & c2x = o.fd_cx.second;
    const BlkCoeffY & c1y = o.fd_cy.first, & c2y = o.fd_cy.second;
    const BlkCoeffZ & c1z = o.fd_cz.first, & c2z = o.fd_cz.second;

    for (int iz=0; iz<FluidBlock::sizeZ; ++iz)
    for (int iy=0; iy<FluidBlock::sizeY; ++iy)
    for (int ix=0; ix<FluidBlock::sizeX; ++ix)
    {
      const FluidElement &L =lab(ix,iy,iz);
      const FluidElement &LW=lab(ix-1,iy,iz), &LE=lab(ix+1,iy,iz);
      const FluidElement &LS=lab(ix,iy-1,iz), &LN=lab(ix,iy+1,iz);
      const FluidElement &LF=lab(ix,iy,iz-1), &LB=lab(ix,iy,iz+1);
      const Real d1udx1 = __FD_2ND(ix, c1x, LW.u, L.u, LE.u);
      const Real d1udy1 = __FD_2ND(iy, c1y, LS.u, L.u, LN.u);
      const Real d1udz1 = __FD_2ND(iz, c1z, LF.u, L.u, LB.u);
      const Real d1vdx1 = __FD_2ND(ix, c1x, LW.v, L.v, LE.v);
      const Real d1vdy1 = __FD_2ND(iy, c1y, LS.v, L.v, LN.v);
      const Real d1vdz1 = __FD_2ND(iz, c1z, LF.v, L.v, LB.v);
      const Real d1wdx1 = __FD_2ND(ix, c1x, LW.w, L.w, LE.w);
      const Real d1wdy1 = __FD_2ND(iy, c1y, LS.w, L.w, LN.w);
      const Real d1wdz1 = __FD_2ND(iz, c1z, LF.w, L.w, LB.w);
      const Real d2udx2 = __FD_2ND(ix, c2x, LW.u, L.u, LE.u);
      const Real d2udy2 = __FD_2ND(iy, c2y, LS.u, L.u, LN.u);
      const Real d2udz2 = __FD_2ND(iz, c2z, LF.u, L.u, LB.u);
      const Real d2vdx2 = __FD_2ND(ix, c2x, LW.v, L.v, LE.v);
      const Real d2vdy2 = __FD_2ND(iy, c2y, LS.v, L.v, LN.v);
      const Real d2vdz2 = __FD_2ND(iz, c2z, LF.v, L.v, LB.v);
      const Real d2wdx2 = __FD_2ND(ix, c2x, LW.w, L.w, LE.w);
      const Real d2wdy2 = __FD_2ND(iy, c2y, LS.w, L.w, LN.w);
      const Real d2wdz2 = __FD_2ND(iz, c2z, LF.w, L.w, LB.w);
      const Real u = L.u+uInf[0], v = L.v+uInf[1], w = L.w+uInf[2];
      const Real duD = d2udx2 + d2udy2 + d2udz2;
      const Real dvD = d2vdx2 + d2vdy2 + d2vdz2;
      const Real dwD = d2wdx2 + d2wdy2 + d2wdz2;
      const Real duA = u * d1udx1 + v * d1udy1 + w * d1udz1;
      const Real dvA = u * d1vdx1 + v * d1vdy1 + w * d1vdz1;
      const Real dwA = u * d1wdx1 + v * d1wdy1 + w * d1wdz1;
      o(ix,iy,iz).tmpU = L.u - dt*duA + dt*mu*duD;
      o(ix,iy,iz).tmpV = L.v - dt*dvA + dt*mu*dvD;
      o(ix,iy,iz).tmpW = L.w - dt*dwA + dt*mu*dwD;
    }
  }
};

void AdvectionDiffusion::operator()(const double dt)
{
  sim.startProfiler("AdvDiff Kernel");
  {
    if(sim.bUseStretchedGrid)
    {
      const KernelAdvectDiffuse_nonUniform K(dt, sim.nu, sim.uinf, sim.lambda);
      compute<KernelAdvectDiffuse_nonUniform>(K);
    }
    else
    {
      //using K_t = KernelAdvectDiffuse_HighOrder;
      using K_t = KernelAdvectDiffuse;
      const K_t K(dt, sim.nu, sim.uinf, sim.lambda);
      compute<K_t>(K);
    }
  }
  sim.stopProfiler();

  sim.startProfiler("AdvDiff copy");
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
  sim.stopProfiler();
  check("AdvectionDiffusion");
}

CubismUP_3D_NAMESPACE_END
