//
//  CubismUP_3D
//  Copyright (c) 2020 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Michalis Chatzimanolakis (michaich@ethz.ch).
//

#pragma once

#include "../SimulationData.h"
#include "../operators/Operator.h"
#include "PoissonSolverBase.h"
#include "PoissonSolverAMRKernels.h"
#include <Cubism/BlockInfo.h>
#include <vector>
#include <cassert>
#include <cstring>

namespace cubismup3d {

static void getZImplParallel(SimulationData &sim, const std::vector<cubism::BlockInfo>& vInfo)
{
  using namespace cubismup3d::poisson_kernels;
  const size_t Nblocks = vInfo.size();

  // We could enable this, we don't really care about denormals.
  // _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);

  // A struct to enforce relative alignment between matrices. The relative
  // alignment of Ax and r MUST NOT be a multiple of 4KB due to cache bank
  // conflicts. See "Haswell and Broadwell pipeline", section
  // "Cache and memory access" here:
  // https://www.agner.org/optimize/microarchitecture.pdf
  struct Tmp {
    // It seems like some offsets with respect to the page boundary of 4KB are
    // faster than the others. (This is accomplished by adding an offset here
    // and using alignas(4096) below). However, this is likely CPU-dependent,
    // so we don't hardcode such fine-tunings here.
    // char offset[0xec0];
    Block r;
    // Ensure p[0+1][0+1][0+xPad] is 64B-aligned for AVX-512 to work.
    char padding1[64 - xPad * sizeof(Real)];
    PaddedBlock p;
    char padding2[xPad * sizeof(Real)];
    Block Ax;
  };
  alignas(64) Tmp tmp{};  // See the kernels cpp file for required alignments.
  Block &r = tmp.r;
  Block &Ax = tmp.Ax;
  PaddedBlock &p = tmp.p;

  #pragma omp for
  for (size_t i = 0; i < Nblocks; ++i) {
    static_assert(sizeof(ScalarBlock) == sizeof(Block));
    assert((uintptr_t)vInfo[i].ptrBlock % kBlockAlignment == 0);
    Block &block = *(Block *)__builtin_assume_aligned(
        vInfo[i].ptrBlock, kBlockAlignment);

    const Real invh = 1 / vInfo[i].h;
    Real rrPartial[NX] = {};
    for (int iz = 0; iz < NZ; ++iz)
    for (int iy = 0; iy < NY; ++iy)
    for (int ix = 0; ix < NX; ++ix) {
      r[iz][iy][ix] = invh * block[iz][iy][ix];
      rrPartial[ix] += r[iz][iy][ix] * r[iz][iy][ix];
      p[iz + 1][iy + 1][ix + xPad] = r[iz][iy][ix];
      block[iz][iy][ix] = 0;
    }
    Real rr = sum(rrPartial);

    const Real sqrNorm0 = (Real)1 / (N * N) * rr;

    if (sqrNorm0 < 1e-32)
      continue;

    const Real *pW = &p[0][0][0] - 1;
    const Real *pE = &p[0][0][0] + 1;

    for (int k = 0; k < 100; ++k) {
      // rr = kernelPoissonGetZInnerReference(p,Ax, r, block, sqrNorm0, rr);
      rr = kernelPoissonGetZInner(p, pW, pE, Ax, r, block, sqrNorm0, rr);
      if (rr <= 0)
        break;
    }
  }
}

class ComputeLHS : public Operator
{
  struct KernelLHSPoisson
  {
    const SimulationData & sim;
    KernelLHSPoisson(const SimulationData&s) : sim(s) {}
    const std::vector<BlockInfo> & lhsInfo = sim.lhsInfo();
    const int Nx = ScalarBlock::sizeX;
    const int Ny = ScalarBlock::sizeY;
    const int Nz = ScalarBlock::sizeZ;
    const StencilInfo stencil{-1,-1,-1,2,2,2,false,{0}};

    void operator()(const ScalarLab & lab, const BlockInfo& info) const
    {
      ScalarBlock & __restrict__ o  = (*sim.lhs)(info.blockID);
      const Real h = info.h; 
      for(int z=0; z<Nz; ++z)
      for(int y=0; y<Ny; ++y)
      for(int x=0; x<Nx; ++x)
      {
        o(x,y,z)   = h*( lab(x-1,y,z) + lab(x+1,y,z) + 
                         lab(x,y-1,z) + lab(x,y+1,z) +
                         lab(x,y,z-1) + lab(x,y,z+1) - 6.0*lab(x,y,z));
      }

      BlockCase<ScalarBlock> * tempCase = (BlockCase<ScalarBlock> *)(lhsInfo[info.blockID].auxiliary);

      if (tempCase == nullptr) return; //no flux corrections needed for this block

      ScalarElement * const faceXm = tempCase -> storedFace[0] ?  & tempCase -> m_pData[0][0] : nullptr;
      ScalarElement * const faceXp = tempCase -> storedFace[1] ?  & tempCase -> m_pData[1][0] : nullptr;
      ScalarElement * const faceYm = tempCase -> storedFace[2] ?  & tempCase -> m_pData[2][0] : nullptr;
      ScalarElement * const faceYp = tempCase -> storedFace[3] ?  & tempCase -> m_pData[3][0] : nullptr;
      ScalarElement * const faceZm = tempCase -> storedFace[4] ?  & tempCase -> m_pData[4][0] : nullptr;
      ScalarElement * const faceZp = tempCase -> storedFace[5] ?  & tempCase -> m_pData[5][0] : nullptr;

      if (faceXm != nullptr)
      {
        const int x = 0;
        for(int z=0; z<Nz; ++z)
        for(int y=0; y<Ny; ++y)
          faceXm[y + Ny*z] = h*(lab(x,y,z) - lab(x-1,y,z));
      }
      if (faceXp != nullptr)
      {
        const int x = Nx-1;
        for(int z=0; z<Nz; ++z)
        for(int y=0; y<Ny; ++y)
          faceXp[y + Ny*z] = h*(lab(x,y,z) - lab(x+1,y,z));
      }
      if (faceYm != nullptr)
      {
        const int y = 0;
        for(int z=0; z<Nz; ++z)
        for(int x=0; x<Nx; ++x)
          faceYm[x + Nx*z] = h*(lab(x,y,z) - lab(x,y-1,z));
      }
      if (faceYp != nullptr)
      {
        const int y = Ny-1;
        for(int z=0; z<Nz; ++z)
        for(int x=0; x<Nx; ++x)
          faceYp[x + Nx*z] = h*(lab(x,y,z) - lab(x,y+1,z));
      }
      if (faceZm != nullptr)
      {
        const int z = 0;
        for(int y=0; y<Ny; ++y)
        for(int x=0; x<Nx; ++x)
          faceZm[x + Nx*y] = h*(lab(x,y,z) - lab(x,y,z-1));
      }
      if (faceZp != nullptr)
      {
        const int z = Nz-1;
        for(int y=0; y<Ny; ++y)
        for(int x=0; x<Nx; ++x)
          faceZp[x + Nx*y] = h*(lab(x,y,z) - lab(x,y,z+1));
      }

    }
  };
  public:
  ComputeLHS(SimulationData & s) : Operator(s) { }
  void operator()(const Real dt)
  {
    compute<ScalarLab>(KernelLHSPoisson(sim),sim.pres,sim.lhs);

    const std::vector<BlockInfo> & vInfo_lhs = sim.lhsInfo();
    const std::vector<BlockInfo> & vInfo_z   = sim.presInfo();
    const int Nx = ScalarBlock::sizeX;
    const int Ny = ScalarBlock::sizeY;
    const int Nz = ScalarBlock::sizeZ;

    if (sim.bMeanConstraint == 0) return;

    if (sim.bMeanConstraint <= 2)
    {
       Real avgP = 0;
       int index = -1;
       #pragma omp parallel for reduction(+ : avgP)
       for(size_t i=0; i<vInfo_lhs.size(); ++i)
       {
          const ScalarBlock & __restrict__ Z  = (*sim.pres)(i);
          const Real h3 = vInfo_z[i].h*vInfo_z[i].h*vInfo_z[i].h;
          if (vInfo_z[i].index[0] == 0 && vInfo_z[i].index[1] == 0 &&  vInfo_z[i].index[2] == 0) index = i;
          for(int z=0; z<Nz; ++z)
          for(int y=0; y<Ny; ++y)
          for(int x=0; x<Nx; ++x)
            avgP += Z(x,y,z).s*h3;
      }
      MPI_Allreduce(MPI_IN_PLACE, &avgP, 1, MPI_Real, MPI_SUM, sim.comm);

      if (sim.bMeanConstraint == 1 && index != -1)
      {
         ScalarBlock & __restrict__ LHS  = (*sim.lhs)(index);
         LHS(0,0,0).s = avgP;
      }
      else if (sim.bMeanConstraint == 2)
      {
         #pragma omp parallel for
         for(size_t i=0; i<vInfo_lhs.size(); ++i)
	 {
            ScalarBlock & __restrict__ LHS = (*sim.lhs)(i);
            const Real h3 = vInfo_lhs[i].h*vInfo_lhs[i].h*vInfo_lhs[i].h;
            for(int z=0; z<Nz; ++z)
            for(int y=0; y<Ny; ++y)
            for(int x=0; x<Nx; ++x)
               LHS(x,y,z).s += avgP*h3;
	 }
      }
    }
    else // > 2
    {
       #pragma omp parallel for
       for(size_t i=0; i<vInfo_lhs.size(); ++i)
       {
          ScalarBlock & __restrict__ LHS = (*sim.lhs)(i);
          const ScalarBlock & __restrict__ Z = (*sim.pres)(i);
          if (vInfo_lhs[i].index[0] == 0 && vInfo_lhs[i].index[1] == 0 && vInfo_lhs[i].index[2] == 0) LHS(0,0,0).s = Z(0,0,0).s;
      }
    }
  }
  std::string getName() { return "ComputeLHS"; }
};

class PoissonSolverAMR : public PoissonSolverBase
{
 protected:
  SimulationData & sim;
  ComputeLHS findLHS;
  void getZ()
  {
    #pragma omp parallel
    {
      getZImplParallel(sim, sim.presInfo());
    }
  }

  size_t _dest(const BlockInfo &info , const int z, const int y, const int x) const
  {
    return BlockType::sizeX * ( BlockType::sizeY * (info.blockID * BlockType::sizeZ  + z) + y) + x;
  }
 public:
  PoissonSolverAMR(SimulationData& s): sim(s),findLHS(s){}
  PoissonSolverAMR(const PoissonSolverAMR& c) = delete; 
  void solve();
};

}//namespace cubismup3d
