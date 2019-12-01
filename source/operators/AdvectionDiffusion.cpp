//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#include "AdvectionDiffusion.h"

CubismUP_3D_NAMESPACE_BEGIN
using namespace cubism;

static constexpr Real EPS = std::numeric_limits<Real>::epsilon();

namespace {

struct KernelAdvectDiffuseBase
{
  KernelAdvectDiffuseBase(const SimulationData&s, double _dt): sim(s),dt(_dt) {}

  const SimulationData & sim;
  const Real dt, mu = sim.nu;
  const std::array<Real, 3>& uInf = sim.uinf;
  const Real fac = std::min((Real)1, sim.uMax_measured * sim.dt / sim.hmean);
  const Real norUinf = std::max({std::fabs(uInf[0]), std::fabs(uInf[1]),
                                 std::fabs(uInf[2]), EPS});
  const Real fadeW = 1 - fac * std::pow(std::max(uInf[0],(Real)0) / norUinf, 2);
  const Real fadeS = 1 - fac * std::pow(std::max(uInf[1],(Real)0) / norUinf, 2);
  const Real fadeF = 1 - fac * std::pow(std::max(uInf[2],(Real)0) / norUinf, 2);
  const Real fadeE = 1 - fac * std::pow(std::min(uInf[0],(Real)0) / norUinf, 2);
  const Real fadeN = 1 - fac * std::pow(std::min(uInf[1],(Real)0) / norUinf, 2);
  const Real fadeB = 1 - fac * std::pow(std::min(uInf[2],(Real)0) / norUinf, 2);
  static constexpr int BEG = -1, END = CUP_BLOCK_SIZE;
  //static constexpr std::array<int, 3> stencil_start = {-1,-1,-1};
  //static constexpr std::array<int, 3> stencil_end   = { 2, 2, 2};
  const StencilInfo stencil{-1,-1,-1, 2,2,2, false, {FE_U,FE_V,FE_W}};

  void applyBCwest(const BlockInfo & I, Lab & L) const {
    if (sim.BCx_flag == wall || sim.BCx_flag == periodic) return;
    else if (I.index[0] not_eq 0) return; // not near boundary
    else if (fadeW >= 1) return; // no momentum killing at this boundary
    assert(fadeW <= 1 && fadeW >= 0);
    for (int iz=-1; iz<=FluidBlock::sizeZ; ++iz)
    for (int iy=-1; iy<=FluidBlock::sizeY; ++iy) {
      L(BEG,iy,iz).u *= fadeW; L(BEG,iy,iz).v *= fadeW; L(BEG,iy,iz).w *= fadeW;
    }
  }

  void applyBCeast(const BlockInfo & I, Lab & L) const {
    if (sim.BCx_flag == wall || sim.BCx_flag == periodic) return;
    else if (I.index[0] not_eq sim.bpdx - 1) return; // not near boundary
    else if (fadeE >= 1) return; // no momentum killing at this boundary
    assert(fadeE <= 1 && fadeE >= 0);
    for (int iz=-1; iz<=FluidBlock::sizeZ; ++iz)
    for (int iy=-1; iy<=FluidBlock::sizeY; ++iy) {
      L(END,iy,iz).u *= fadeE; L(END,iy,iz).v *= fadeE; L(END,iy,iz).w *= fadeE;
    }
  }

  void applyBCsouth(const BlockInfo & I, Lab & L) const {
    if (sim.BCy_flag == wall || sim.BCy_flag == periodic) return;
    else if (I.index[1] not_eq 0) return; // not near boundary
    else if (fadeS >= 1) return; // no momentum killing at this boundary
    assert(fadeS <= 1 && fadeS >= 0);
    for (int iz=-1; iz<=FluidBlock::sizeZ; ++iz)
    for (int ix=-1; ix<=FluidBlock::sizeX; ++ix) {
      L(ix,BEG,iz).u *= fadeS; L(ix,BEG,iz).v *= fadeS; L(ix,BEG,iz).w *= fadeS;
    }
  }

  void applyBCnorth(const BlockInfo & I, Lab & L) const {
    if (sim.BCy_flag == wall || sim.BCy_flag == periodic) return;
    else if (I.index[1] not_eq sim.bpdy - 1) return; // not near boundary
    else if (fadeN >= 1) return; // no momentum killing at this boundary
    assert(fadeN <= 1 && fadeN >= 0);
    for (int iz=-1; iz<=FluidBlock::sizeZ; ++iz)
    for (int ix=-1; ix<=FluidBlock::sizeX; ++ix) {
      L(ix,END,iz).u *= fadeN; L(ix,END,iz).v *= fadeN; L(ix,END,iz).w *= fadeN;
    }
  }

  void applyBCfront(const BlockInfo & I, Lab & L) const {
    if (sim.BCz_flag == wall || sim.BCz_flag == periodic) return;
    else if (I.index[2] not_eq 0) return; // not near boundary
    else if (fadeF >= 1) return; // no momentum killing at this boundary
    assert(fadeF <= 1 && fadeF >= 0);
    for (int iy=-1; iy<=FluidBlock::sizeY; ++iy)
    for (int ix=-1; ix<=FluidBlock::sizeX; ++ix) {
      L(ix,iy,BEG).u *= fadeF; L(ix,iy,BEG).v *= fadeF; L(ix,iy,BEG).w *= fadeF;
    }
  }

  void applyBCback(const BlockInfo & I, Lab & L) const {
    if (sim.BCz_flag == wall || sim.BCz_flag == periodic) return;
    else if (I.index[2] not_eq sim.bpdz - 1) return; // not near boundary
    else if (fadeB >= 1) return; // no momentum killing at this boundary
    assert(fadeB <= 1 && fadeB >= 0);
    for (int iy=-1; iy<=FluidBlock::sizeY; ++iy)
    for (int ix=-1; ix<=FluidBlock::sizeX; ++ix) {
      L(ix,iy,END).u *= fadeB; L(ix,iy,END).v *= fadeB; L(ix,iy,END).w *= fadeB;
    }
  }
};

struct KernelAdvectDiffuse : public KernelAdvectDiffuseBase
{
  KernelAdvectDiffuse(const SimulationData&s, double _dt) :
    KernelAdvectDiffuseBase(s, _dt) {}

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
    KernelAdvectDiffuseBase(s, _dt) {}

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
    const KernelAdvectDiffuse K(sim, dt);
    compute(K);
    sim.stopProfiler();
    sim.startProfiler("AdvDiff copy");
    //const UpdateAndCorrectInflow U(sim);
    const UpdateAndCorrectInflow_nonUniform U(sim);
    U.operate();
    sim.stopProfiler();
  }
  check("AdvectionDiffusion");
}

CubismUP_3D_NAMESPACE_END
