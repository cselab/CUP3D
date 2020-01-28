//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#include "AdvectionDiffusion.h"
#include "../obstacles/ObstacleVector.h"

CubismUP_3D_NAMESPACE_BEGIN
using namespace cubism;

static constexpr Real EPS = std::numeric_limits<Real>::epsilon();

namespace {

struct KernelAdvectDiffuseBase
{
  KernelAdvectDiffuseBase(const SimulationData&s, double _dt, int beg, int end)
  : sim(s), dt(_dt), stencilBeg(beg), stencilEnd(end) {
    //printf("%d %d %e %e %e %e %e %e %e %e\n", loopBeg, loopEnd, CFL,
    //  norUinf, fadeW, fadeS, fadeF, fadeE, fadeN, fadeB);
  }

  const SimulationData & sim;
  const Real dt, mu = sim.nu;
  const int stencilBeg, stencilEnd;
  const std::array<Real, 3>& uInf = sim.uinf;
  const int loopBeg = stencilBeg, loopEnd = CUP_BLOCK_SIZE-1 + stencilEnd;

  const Real CFL = std::min((Real)1, sim.uMax_measured * sim.dt / sim.hmean);
  const Real norUinf = std::max({std::fabs(uInf[0]), std::fabs(uInf[1]),
                                 std::fabs(uInf[2]), EPS});
  const Real fadeW = 1 - CFL * std::pow(std::max(uInf[0],(Real)0) / norUinf, 2);
  const Real fadeS = 1 - CFL * std::pow(std::max(uInf[1],(Real)0) / norUinf, 2);
  const Real fadeF = 1 - CFL * std::pow(std::max(uInf[2],(Real)0) / norUinf, 2);
  const Real fadeE = 1 - CFL * std::pow(std::min(uInf[0],(Real)0) / norUinf, 2);
  const Real fadeN = 1 - CFL * std::pow(std::min(uInf[1],(Real)0) / norUinf, 2);
  const Real fadeB = 1 - CFL * std::pow(std::min(uInf[2],(Real)0) / norUinf, 2);
  const StencilInfo stencil{stencilBeg, stencilBeg, stencilBeg,
                            stencilEnd, stencilEnd, stencilEnd,
                            false, {FE_U,FE_V,FE_W}};

  void applyBCwest(const BlockInfo & I, Lab & L) const
  {
    if (sim.BCx_flag == wall || sim.BCx_flag == periodic) return;
    else if (I.index[0] not_eq 0) return; // not near boundary
    else if (fadeW >= 1) return; // no momentum killing at this boundary
    for (int ix = loopBeg; ix < 0; ++ix)
    {
      const Real fac = std::pow(fadeW, 0 - ix);
      assert(fac <= 1 && fac >= 0);
      for (int iz = loopBeg; iz < loopEnd; ++iz)
      for (int iy = loopBeg; iy < loopEnd; ++iy)
      {
        L(ix,iy,iz).u *= fac; L(ix,iy,iz).v *= fac; L(ix,iy,iz).w *= fac;
      }
    }
  }

  void applyBCeast(const BlockInfo & I, Lab & L) const
  {
    if (sim.BCx_flag == wall || sim.BCx_flag == periodic) return;
    else if (I.index[0] not_eq sim.bpdx - 1) return; // not near boundary
    else if (fadeE >= 1) return; // no momentum killing at this boundary
    for (int ix = CUP_BLOCK_SIZE; ix < loopEnd; ++ix)
    {
      const Real fac = std::pow(fadeE, ix - CUP_BLOCK_SIZE + 1);
      assert(fac <= 1 && fac >= 0);
      for (int iz = loopBeg; iz < loopEnd; ++iz)
      for (int iy = loopBeg; iy < loopEnd; ++iy)
      {
        L(ix,iy,iz).u *= fac; L(ix,iy,iz).v *= fac; L(ix,iy,iz).w *= fac;
      }
    }
  }

  void applyBCsouth(const BlockInfo & I, Lab & L) const
  {
    if (sim.BCy_flag == wall || sim.BCy_flag == periodic) return;
    else if (I.index[1] not_eq 0) return; // not near boundary
    else if (fadeS >= 1) return; // no momentum killing at this boundary
    for (int iy = loopBeg; iy < 0; ++iy)
    {
      const Real fac = std::pow(fadeS, 0 - iy);
      assert(fac <= 1 && fac >= 0);
      for (int iz = loopBeg; iz < loopEnd; ++iz)
      for (int ix = loopBeg; ix < loopEnd; ++ix)
      {
        L(ix,iy,iz).u *= fac; L(ix,iy,iz).v *= fac; L(ix,iy,iz).w *= fac;
      }
    }
  }

  void applyBCnorth(const BlockInfo & I, Lab & L) const
  {
    if (sim.BCy_flag == wall || sim.BCy_flag == periodic) return;
    else if (I.index[1] not_eq sim.bpdy - 1) return; // not near boundary
    else if (fadeN >= 1) return; // no momentum killing at this boundary
    for (int iy = CUP_BLOCK_SIZE; iy < loopEnd; ++iy)
    {
      const Real fac = std::pow(fadeN, iy - CUP_BLOCK_SIZE + 1);
      assert(fac <= 1 && fac >= 0);
      for (int iz = loopBeg; iz < loopEnd; ++iz)
      for (int ix = loopBeg; ix < loopEnd; ++ix)
      {
        L(ix,iy,iz).u *= fac; L(ix,iy,iz).v *= fac; L(ix,iy,iz).w *= fac;
      }
    }
  }

  void applyBCfront(const BlockInfo & I, Lab & L) const
  {
    if (sim.BCz_flag == wall || sim.BCz_flag == periodic) return;
    else if (I.index[2] not_eq 0) return; // not near boundary
    else if (fadeF >= 1) return; // no momentum killing at this boundary
    for (int iz = loopBeg; iz < 0; ++iz)
    {
      const Real fac = std::pow(fadeF, 0 - iz);
      assert(fac <= 1 && fac >= 0);
      for (int iy = loopBeg; iy < loopEnd; ++iy)
      for (int ix = loopBeg; ix < loopEnd; ++ix)
      {
        L(ix,iy,iz).u *= fac; L(ix,iy,iz).v *= fac; L(ix,iy,iz).w *= fac;
      }
    }
  }

  void applyBCback(const BlockInfo & I, Lab & L) const
  {
    if (sim.BCz_flag == wall || sim.BCz_flag == periodic) return;
    else if (I.index[2] not_eq sim.bpdz - 1) return; // not near boundary
    else if (fadeB >= 1) return; // no momentum killing at this boundary
    for (int iz = CUP_BLOCK_SIZE; iz < loopEnd; ++iz)
    {
      const Real fac = std::pow(fadeB, iz - CUP_BLOCK_SIZE + 1);
      assert(fac <= 1 && fac >= 0);
      for (int iy = loopBeg; iy < loopEnd; ++iy)
      for (int ix = loopBeg; ix < loopEnd; ++ix)
      {
        L(ix,iy,iz).u *= fac; L(ix,iy,iz).v *= fac; L(ix,iy,iz).w *= fac;
      }
    }
  }
};

struct KernelAdvectDiffuse3rdOrderUpwind : public KernelAdvectDiffuseBase
{
  KernelAdvectDiffuse3rdOrderUpwind(const SimulationData&s, double _dt) :
    KernelAdvectDiffuseBase(s, _dt, -2, 3) {}

  template <typename Lab, typename BlockType>
  void operator()(Lab & lab, const BlockInfo& info, BlockType& o) const
  {
    const Real facA = -dt/(6*info.h_gridpoint);
    const Real facD = (mu/info.h_gridpoint) * (dt/info.h_gridpoint);
    applyBCwest(info, lab);
    applyBCeast(info, lab);
    applyBCsouth(info, lab);
    applyBCnorth(info, lab);
    applyBCfront(info, lab);
    applyBCback(info, lab);

    for (int iz=0; iz<FluidBlock::sizeZ; ++iz)
    for (int iy=0; iy<FluidBlock::sizeY; ++iy)
    for (int ix=0; ix<FluidBlock::sizeX; ++ix)
    {
      const FluidElement &L   = lab(ix,iy,iz);
      const FluidElement &LW  = lab(ix-1,iy,iz), &LE  = lab(ix+1,iy,iz);
      const FluidElement &LS  = lab(ix,iy-1,iz), &LN  = lab(ix,iy+1,iz);
      const FluidElement &LF  = lab(ix,iy,iz-1), &LB  = lab(ix,iy,iz+1);
      const FluidElement &LW2 = lab(ix-2,iy,iz), &LE2 = lab(ix+2,iy,iz);
      const FluidElement &LS2 = lab(ix,iy-2,iz), &LN2 = lab(ix,iy+2,iz);
      const FluidElement &LF2 = lab(ix,iy,iz-2), &LB2 = lab(ix,iy,iz+2);

      const Real u = L.u+uInf[0], v = L.v+uInf[1], w = L.w+uInf[2];
      #if 1
      const Real dudxM =         2*LE.u +3*L.u -6*LW.u +LW2.u;
      const Real dudxP = -LE2.u +6*LE.u -3*L.u -2*LW.u;
      const Real dvdxM =         2*LE.v +3*L.v -6*LW.v +LW2.v;
      const Real dvdxP = -LE2.v +6*LE.v -3*L.v -2*LW.v;
      const Real dwdxM =         2*LE.w +3*L.w -6*LW.w +LW2.w;
      const Real dwdxP = -LE2.w +6*LE.w -3*L.w -2*LW.w;
      const Real dudyM =         2*LN.u +3*L.u -6*LS.u +LS2.u;
      const Real dudyP = -LN2.u +6*LN.u -3*L.u -2*LS.u;
      const Real dvdyM =         2*LN.v +3*L.v -6*LS.v +LS2.v;
      const Real dvdyP = -LN2.v +6*LN.v -3*L.v -2*LS.v;
      const Real dwdyM =         2*LN.w +3*L.w -6*LS.w +LS2.w;
      const Real dwdyP = -LN2.w +6*LN.w -3*L.w -2*LS.w;
      const Real dudzM =         2*LB.u +3*L.u -6*LF.u +LF2.u;
      const Real dudzP = -LB2.u +6*LB.u -3*L.u -2*LF.u;
      const Real dvdzM =         2*LB.v +3*L.v -6*LF.v +LF2.v;
      const Real dvdzP = -LB2.v +6*LB.v -3*L.v -2*LF.v;
      const Real dwdzM =         2*LB.w +3*L.w -6*LF.w +LF2.w;
      const Real dwdzP = -LB2.w +6*LB.w -3*L.w -2*LF.w;
      const Real dudxC = LE.u-LW.u, dvdxC = LE.v-LW.v, dwdxC = LE.w-LW.w;
      const Real dudyC = LN.u-LS.u, dvdyC = LN.v-LS.v, dwdyC = LN.w-LS.w;
      const Real dudzC = LB.u-LF.u, dvdzC = LB.v-LF.v, dwdzC = LB.w-LF.w;
      #if 0
      const Real invU = 1.0 / std::sqrt(u*u + v*v + w*w);
      const Real UP = std::max((Real)0, u), UM = std::min((Real)0, u);
      const Real VP = std::max((Real)0, v), VM = std::min((Real)0, v);
      const Real WP = std::max((Real)0, w), WM = std::min((Real)0, w);
      const Real WXM = UP*invU, WXP = -UM*invU, WXC = 1-std::fabs(u)*invU;
      const Real WYM = VP*invU, WYP = -VM*invU, WYC = 1-std::fabs(v)*invU;
      const Real WZM = WP*invU, WZP = -WM*invU, WZC = 1-std::fabs(w)*invU;
      #else
      const Real invU = 1.0 / (u*u + v*v + w*w);
      const Real UP = std::max((Real)0, u), UM = std::min((Real)0, u);
      const Real VP = std::max((Real)0, v), VM = std::min((Real)0, v);
      const Real WP = std::max((Real)0, w), WM = std::min((Real)0, w);
      const Real WXM = UP*UP*invU, WXP = UM*UM*invU, WXC = 1-u*u*invU;
      const Real WYM = VP*VP*invU, WYP = VM*VM*invU, WYC = 1-v*v*invU;
      const Real WZM = WP*WP*invU, WZP = WM*WM*invU, WZC = 1-w*w*invU;
      #endif
      //printf("%e %e %e\n", WXM+WXP+WXC-1,WYM+WYP+WYC-1,WZM+WZP+WZC-1);
      const Real dudx = WXM * dudxM + WXP * dudxP + WXC * 3*dudxC;
      const Real dvdx = WXM * dvdxM + WXP * dvdxP + WXC * 3*dvdxC;
      const Real dwdx = WXM * dwdxM + WXP * dwdxP + WXC * 3*dwdxC;
      const Real dudy = WYM * dudyM + WYP * dudyP + WYC * 3*dudyC;
      const Real dvdy = WYM * dvdyM + WYP * dvdyP + WYC * 3*dvdyC;
      const Real dwdy = WYM * dwdyM + WYP * dwdyP + WYC * 3*dwdyC;
      const Real dudz = WZM * dudzM + WZP * dudzP + WZC * 3*dudzC;
      const Real dvdz = WZM * dvdzM + WZP * dvdzP + WZC * 3*dvdzC;
      const Real dwdz = WZM * dwdzM + WZP * dwdzP + WZC * 3*dwdzC;
      #else
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

struct KernelAdvectDiffuse : public KernelAdvectDiffuseBase
{
  KernelAdvectDiffuse(const SimulationData&s, double _dt) :
    KernelAdvectDiffuseBase(s, _dt, -1, 2) {}

  template <typename Lab, typename BlockType>
  void operator()(Lab & lab, const BlockInfo& info, BlockType& o) const
  {
    const Real facA = -dt/(2*info.h_gridpoint);
    const Real facD = (mu/info.h_gridpoint) * (dt/info.h_gridpoint);
    applyBCwest(info, lab);
    applyBCeast(info, lab);
    applyBCsouth(info, lab);
    applyBCnorth(info, lab);
    applyBCfront(info, lab);
    applyBCback(info, lab);

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
      o(ix,iy,iz).tmpU = L.u + facA*duA + facD*duD;
      o(ix,iy,iz).tmpV = L.v + facA*dvA + facD*dvD;
      o(ix,iy,iz).tmpW = L.w + facA*dwA + facD*dwD;
    }
  }
};

struct KernelAdvectDiffuse_nonUniform : public KernelAdvectDiffuseBase
{
  KernelAdvectDiffuse_nonUniform(const SimulationData&s, double _dt) :
    KernelAdvectDiffuseBase(s, _dt, -1, 2) {}

  template <typename Lab, typename BlockType>
  void operator()(Lab & lab, const BlockInfo& info, BlockType& o) const
  {
    // FD coefficients for first and second derivative
    const BlkCoeffX & c1x = o.fd_cx.first, & c2x = o.fd_cx.second;
    const BlkCoeffY & c1y = o.fd_cy.first, & c2y = o.fd_cy.second;
    const BlkCoeffZ & c1z = o.fd_cz.first, & c2z = o.fd_cz.second;
    applyBCwest(info, lab);
    applyBCeast(info, lab);
    applyBCsouth(info, lab);
    applyBCnorth(info, lab);
    applyBCfront(info, lab);
    applyBCback(info, lab);

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

struct UpdateAndCorrectInflow
{
  SimulationData & sim;
  FluidGridMPI * const grid = sim.grid;
  const std::vector<cubism::BlockInfo>& vInfo = sim.vInfo();

  static constexpr int BEG = 0, END = CUP_BLOCK_SIZE-1;
  inline bool isW(const BlockInfo&I) const {
    if (sim.BCx_flag == wall || sim.BCx_flag == periodic) return false;
    return I.index[0] == 0;
  };
  inline bool isE(const BlockInfo&I) const {
    if (sim.BCx_flag == wall || sim.BCx_flag == periodic) return false;
    return I.index[0] == sim.bpdx-1;
  };
  inline bool isS(const BlockInfo&I) const {
    if (sim.BCy_flag == wall || sim.BCy_flag == periodic) return false;
    return I.index[1] == 0;
  };
  inline bool isN(const BlockInfo&I) const {
    if (sim.BCy_flag == wall || sim.BCy_flag == periodic) return false;
    return I.index[1] == sim.bpdy-1;
  };
  inline bool isF(const BlockInfo&I) const {
    if (sim.BCz_flag == wall || sim.BCz_flag == periodic) return false;
    return I.index[2] == 0;
  };
  inline bool isB(const BlockInfo&I) const {
    if (sim.BCz_flag == wall || sim.BCz_flag == periodic) return false;
    return I.index[2] == sim.bpdz-1;
  };

  UpdateAndCorrectInflow(SimulationData & s) : sim(s) { }
  virtual ~UpdateAndCorrectInflow() {}

  virtual void operate() const
  {
    double sumInflow = 0;
    #pragma omp parallel for schedule(static) reduction(+:sumInflow)
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

      if (isW(vInfo[i])) for (int iz=0; iz<FluidBlock::sizeZ; ++iz)
        for (int iy=0; iy<FluidBlock::sizeY; ++iy) sumInflow -= b(BEG,iy,iz).u;

      if (isE(vInfo[i])) for (int iz=0; iz<FluidBlock::sizeZ; ++iz)
        for (int iy=0; iy<FluidBlock::sizeY; ++iy) sumInflow += b(END,iy,iz).u;

      if (isS(vInfo[i])) for (int iz=0; iz<FluidBlock::sizeZ; ++iz)
        for (int ix=0; ix<FluidBlock::sizeX; ++ix) sumInflow -= b(ix,BEG,iz).v;

      if (isN(vInfo[i])) for (int iz=0; iz<FluidBlock::sizeZ; ++iz)
        for (int ix=0; ix<FluidBlock::sizeX; ++ix) sumInflow += b(ix,END,iz).v;

      if (isF(vInfo[i])) for (int iy=0; iy<FluidBlock::sizeY; ++iy)
        for (int ix=0; ix<FluidBlock::sizeX; ++ix) sumInflow -= b(ix,iy,BEG).w;

      if (isB(vInfo[i])) for (int iy=0; iy<FluidBlock::sizeY; ++iy)
        for (int ix=0; ix<FluidBlock::sizeX; ++ix) sumInflow += b(ix,iy,END).w;
    }

    const auto & comm = grid->getCartComm();
    MPI_Allreduce(MPI_IN_PLACE, & sumInflow, 1, MPI_DOUBLE, MPI_SUM, comm);
    const auto nTotX = FluidBlock::sizeX * sim.bpdx;
    const auto nTotY = FluidBlock::sizeY * sim.bpdy;
    const auto nTotZ = FluidBlock::sizeZ * sim.bpdz;
    const Real corr = sumInflow / (2*(nTotX*nTotY + nTotX*nTotZ + nTotY*nTotZ));

    if(std::fabs(corr) < EPS) return;
    if(sim.verbose) printf("Inflow correction %e\n", corr);

    #pragma omp parallel for schedule(static)
    for(size_t i=0; i<vInfo.size(); i++)
    {
      FluidBlock& b = *(FluidBlock*) vInfo[i].ptrBlock;

      if(isW(vInfo[i])) for (int iz=0; iz<FluidBlock::sizeZ; ++iz)
        for (int iy=0; iy<FluidBlock::sizeY; ++iy) b(BEG,iy,iz).u += corr;

      if(isE(vInfo[i])) for (int iz=0; iz<FluidBlock::sizeZ; ++iz)
        for (int iy=0; iy<FluidBlock::sizeY; ++iy) b(END,iy,iz).u -= corr;

      if(isS(vInfo[i])) for (int iz=0; iz<FluidBlock::sizeZ; ++iz)
        for (int ix=0; ix<FluidBlock::sizeX; ++ix) b(ix,BEG,iz).v += corr;

      if(isN(vInfo[i])) for (int iz=0; iz<FluidBlock::sizeZ; ++iz)
        for (int ix=0; ix<FluidBlock::sizeX; ++ix) b(ix,END,iz).v -= corr;

      if(isF(vInfo[i])) for (int iy=0; iy<FluidBlock::sizeY; ++iy)
        for (int ix=0; ix<FluidBlock::sizeX; ++ix) b(ix,iy,BEG).w += corr;

      if(isB(vInfo[i])) for (int iy=0; iy<FluidBlock::sizeY; ++iy)
        for (int ix=0; ix<FluidBlock::sizeX; ++ix) b(ix,iy,END).w -= corr;
    }
  }
};

struct UpdateAndCorrectInflow_nonUniform : public UpdateAndCorrectInflow
{
  UpdateAndCorrectInflow_nonUniform(SimulationData & s) :
    UpdateAndCorrectInflow(s) { }

  void operate() const override
  {
    double sumInflow = 0, throughFlow = 0;
    #pragma omp parallel for schedule(static) reduction(+:sumInflow,throughFlow)
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

      if(isW(vInfo[i])) for (int iz=0; iz<FluidBlock::sizeZ; ++iz)
                        for (int iy=0; iy<FluidBlock::sizeY; ++iy) {
        sumInflow -= b(BEG,iy,iz).u; throughFlow += std::fabs(b(BEG,iy,iz).u);
      }

      if(isE(vInfo[i])) for (int iz=0; iz<FluidBlock::sizeZ; ++iz)
                        for (int iy=0; iy<FluidBlock::sizeY; ++iy) {
        sumInflow += b(END,iy,iz).u; throughFlow += std::fabs(b(END,iy,iz).u);
      }

      if(isS(vInfo[i])) for (int iz=0; iz<FluidBlock::sizeZ; ++iz)
                        for (int ix=0; ix<FluidBlock::sizeX; ++ix) {
        sumInflow -= b(ix,BEG,iz).v; throughFlow += std::fabs(b(ix,BEG,iz).v);
      }

      if(isN(vInfo[i])) for (int iz=0; iz<FluidBlock::sizeZ; ++iz)
                        for (int ix=0; ix<FluidBlock::sizeX; ++ix) {
        sumInflow += b(ix,END,iz).v; throughFlow += std::fabs(b(ix,END,iz).v);
      }

      if(isF(vInfo[i])) for (int iy=0; iy<FluidBlock::sizeY; ++iy)
                        for (int ix=0; ix<FluidBlock::sizeX; ++ix) {
        sumInflow -= b(ix,iy,BEG).w; throughFlow += std::fabs(b(ix,iy,BEG).w);
      }

      if(isB(vInfo[i])) for (int iy=0; iy<FluidBlock::sizeY; ++iy)
                        for (int ix=0; ix<FluidBlock::sizeX; ++ix) {
        sumInflow += b(ix,iy,END).w; throughFlow += std::fabs(b(ix,iy,END).w);
      }
    }

    double sums[2] = {sumInflow, throughFlow};
    MPI_Allreduce(MPI_IN_PLACE, sums,2,MPI_DOUBLE,MPI_SUM, grid->getCartComm());
    const Real corr = sums[0] / std::max((double)EPS, sums[1]);

    if(std::fabs(corr) < EPS) return;
    if(sim.verbose) printf("Relative inflow correction %e\n", corr);

    #pragma omp parallel for schedule(static)
    for(size_t i=0; i<vInfo.size(); i++)
    {
      FluidBlock& b = *(FluidBlock*) vInfo[i].ptrBlock;

      if(isW(vInfo[i])) for (int iz=0; iz<FluidBlock::sizeZ; ++iz)
                        for (int iy=0; iy<FluidBlock::sizeY; ++iy)
        b(BEG,iy,iz).u += corr * std::fabs(b(BEG,iy,iz).u);

      if(isE(vInfo[i])) for (int iz=0; iz<FluidBlock::sizeZ; ++iz)
                        for (int iy=0; iy<FluidBlock::sizeY; ++iy)
        b(END,iy,iz).u -= corr * std::fabs(b(END,iy,iz).u);

      if(isS(vInfo[i])) for (int iz=0; iz<FluidBlock::sizeZ; ++iz)
                        for (int ix=0; ix<FluidBlock::sizeX; ++ix)
        b(ix,BEG,iz).v += corr * std::fabs(b(ix,BEG,iz).v);

      if(isN(vInfo[i])) for (int iz=0; iz<FluidBlock::sizeZ; ++iz)
                        for (int ix=0; ix<FluidBlock::sizeX; ++ix)
        b(ix,END,iz).v -= corr * std::fabs(b(ix,END,iz).v);

      if(isF(vInfo[i])) for (int iy=0; iy<FluidBlock::sizeY; ++iy)
                        for (int ix=0; ix<FluidBlock::sizeX; ++ix)
        b(ix,iy,BEG).w += corr * std::fabs(b(ix,iy,BEG).w);

      if(isB(vInfo[i])) for (int iy=0; iy<FluidBlock::sizeY; ++iy)
                        for (int ix=0; ix<FluidBlock::sizeX; ++ix)
        b(ix,iy,END).w -= corr * std::fabs(b(ix,iy,END).w);
    }
  }
};

}

void AdvectionDiffusion::operator()(const double dt)
{
  if(sim.bUseStretchedGrid)
  {
    sim.startProfiler("AdvDiff Kernel");
    const KernelAdvectDiffuse_nonUniform K(sim, dt);
    compute(K);
    sim.stopProfiler();
    sim.startProfiler("AdvDiff copy");
    const UpdateAndCorrectInflow_nonUniform U(sim);
    U.operate();
    sim.stopProfiler();
  }
  else
  {
    sim.startProfiler("AdvDiff Kernel");
    if(sim.obstacle_vector->nObstacles() == 0) {
      const KernelAdvectDiffuse K(sim, dt);
      //const KernelAdvectDiffuse3rdOrderUpwind K(sim, dt);
      compute(K);
    } else {
      const KernelAdvectDiffuse K(sim, dt);
      compute(K);
    }
    sim.stopProfiler();
    sim.startProfiler("AdvDiff copy");
    const UpdateAndCorrectInflow U(sim);
    //const UpdateAndCorrectInflow_nonUniform U(sim);
    U.operate();
    sim.stopProfiler();
  }
  check("AdvectionDiffusion");
}

CubismUP_3D_NAMESPACE_END
