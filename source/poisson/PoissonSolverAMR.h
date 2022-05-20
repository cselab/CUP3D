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
#include <Cubism/BlockInfo.h>
#include <vector>
#include <cassert>
#include <cstring>

namespace cubismup3d {

class ComputeLHS : public Operator
{
  struct KernelLHSPoisson
  {
    const SimulationData & sim;
    KernelLHSPoisson(const SimulationData&s) : sim(s) {}
    std::vector<BlockInfo> & vInfo_lhs = sim.lhs->getBlocksInfo();
  
    const StencilInfo stencil{-1,-1,-1,2,2,2,false,{0}};
    void operator()(ScalarLab & lab, const BlockInfo& info) const
    {
      ScalarBlock & __restrict__ o  = *(ScalarBlock*) vInfo_lhs[info.blockID].ptrBlock;
      const Real h = info.h; 
      for(int iz=0; iz<ScalarBlock::sizeZ; ++iz)
      for(int iy=0; iy<ScalarBlock::sizeY; ++iy)
      for(int ix=0; ix<ScalarBlock::sizeX; ++ix)
      {
        o(ix,iy,iz).s   = h*( lab(ix-1,iy,iz).s + lab(ix+1,iy,iz).s + 
                              lab(ix,iy-1,iz).s + lab(ix,iy+1,iz).s +
                              lab(ix,iy,iz-1).s + lab(ix,iy,iz+1).s - 6.0*lab(ix,iy,iz).s);
      }

      BlockCase<ScalarBlock> * tempCase = (BlockCase<ScalarBlock> *)(vInfo_lhs[info.blockID].auxiliary);
      ScalarBlock::ElementType * faceXm = nullptr;
      ScalarBlock::ElementType * faceXp = nullptr;
      ScalarBlock::ElementType * faceYm = nullptr;
      ScalarBlock::ElementType * faceYp = nullptr;
      ScalarBlock::ElementType * faceZp = nullptr;
      ScalarBlock::ElementType * faceZm = nullptr;
      if (tempCase != nullptr)
      {
        faceXm = tempCase -> storedFace[0] ?  & tempCase -> m_pData[0][0] : nullptr;
        faceXp = tempCase -> storedFace[1] ?  & tempCase -> m_pData[1][0] : nullptr;
        faceYm = tempCase -> storedFace[2] ?  & tempCase -> m_pData[2][0] : nullptr;
        faceYp = tempCase -> storedFace[3] ?  & tempCase -> m_pData[3][0] : nullptr;
        faceZm = tempCase -> storedFace[4] ?  & tempCase -> m_pData[4][0] : nullptr;
        faceZp = tempCase -> storedFace[5] ?  & tempCase -> m_pData[5][0] : nullptr;
      }
      if (faceXm != nullptr)
      {
        int ix = 0;
        for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
        for(int iy=0; iy<FluidBlock::sizeY; ++iy)
        {
          faceXm[iy + FluidBlock::sizeY * iz].clear();
          faceXm[iy + FluidBlock::sizeY * iz].s = h*(lab(ix,iy,iz).s - lab(ix-1,iy,iz).s);
        }
      }
      if (faceXp != nullptr)
      {
        int ix = FluidBlock::sizeX-1;
        for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
        for(int iy=0; iy<FluidBlock::sizeY; ++iy)
        {
          faceXp[iy + FluidBlock::sizeY * iz].clear();
          faceXp[iy + FluidBlock::sizeY * iz].s = h*(lab(ix,iy,iz).s - lab(ix+1,iy,iz).s);
        }
      }
      if (faceYm != nullptr)
      {
        int iy = 0;
        for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
        for(int ix=0; ix<FluidBlock::sizeX; ++ix)
        {
          faceYm[ix + FluidBlock::sizeX * iz].clear();
          faceYm[ix + FluidBlock::sizeX * iz].s = h*(lab(ix,iy,iz).s - lab(ix,iy-1,iz).s);
        }
      }
      if (faceYp != nullptr)
      {
        int iy = FluidBlock::sizeY-1;
        for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
        for(int ix=0; ix<FluidBlock::sizeX; ++ix)
        {
          faceYp[ix + FluidBlock::sizeX * iz].clear();
          faceYp[ix + FluidBlock::sizeX * iz].s = h*(lab(ix,iy,iz).s - lab(ix,iy+1,iz).s);
        }
      }
      if (faceZm != nullptr)
      {
        int iz = 0;
        for(int iy=0; iy<FluidBlock::sizeY; ++iy)
        for(int ix=0; ix<FluidBlock::sizeX; ++ix)
        {
          faceZm[ix + FluidBlock::sizeX * iy].clear();
          faceZm[ix + FluidBlock::sizeX * iy].s = h*(lab(ix,iy,iz).s - lab(ix,iy,iz-1).s);
        }
      }
      if (faceZp != nullptr)
      {
        int iz = FluidBlock::sizeZ-1;
        for(int iy=0; iy<FluidBlock::sizeY; ++iy)
        for(int ix=0; ix<FluidBlock::sizeX; ++ix)
        {
          faceZp[ix + FluidBlock::sizeX * iy].clear();
          faceZp[ix + FluidBlock::sizeX * iy].s = h*(lab(ix,iy,iz).s - lab(ix,iy,iz+1).s);
        }
      }

    }
  };
  public:
  ComputeLHS(SimulationData & s) : Operator(s) { }
  void operator()(const Real dt)
  {
    std::vector<BlockInfo> & vInfo_lhs = sim.lhs->getBlocksInfo();
    std::vector<BlockInfo> & vInfo_z   = sim.z  ->getBlocksInfo();
    KernelLHSPoisson KPoisson(sim);
    cubism::compute<ScalarLab>(KPoisson,sim.z,sim.lhs);

    if (sim.bMeanConstraint == 0) return;

    if (sim.bMeanConstraint <= 2)
    {
       Real avgP = 0;
       int index = -1;
       #pragma omp parallel for reduction(+ : avgP)
       for(size_t i=0; i<vInfo_lhs.size(); ++i)
       {
          ScalarBlock & __restrict__ Z  = *(ScalarBlock*) vInfo_z[i].ptrBlock;
          const Real h3 = vInfo_z[i].h*vInfo_z[i].h*vInfo_z[i].h;
          if (vInfo_z[i].index[0] == 0 && 
              vInfo_z[i].index[1] == 0 && 
              vInfo_z[i].index[2] == 0)
            index = i;
          for(int iz=0; iz<FluidBlock::sizeZ; iz++)
          for(int iy=0; iy<FluidBlock::sizeY; iy++)
          for(int ix=0; ix<FluidBlock::sizeX; ix++)
            avgP += Z(ix,iy,iz).s*h3;
      }
      MPI_Allreduce(MPI_IN_PLACE, &avgP, 1, MPI_Real, MPI_SUM, sim.grid->getCartComm());

      if (sim.bMeanConstraint == 1 && index != -1)
      {
         ScalarBlock & __restrict__ LHS = *(ScalarBlock*) vInfo_lhs[index].ptrBlock;
         LHS(0,0,0).s = avgP;
      }
      else if (sim.bMeanConstraint == 2)
      {
         #pragma omp parallel for reduction(+ : avgP)
         for(size_t i=0; i<vInfo_lhs.size(); ++i)
	 {
            ScalarBlock & __restrict__ LHS = *(ScalarBlock*) vInfo_lhs[i].ptrBlock;
            const Real h3 = vInfo_lhs[i].h*vInfo_lhs[i].h*vInfo_lhs[i].h;
            for(int iz=0; iz<FluidBlock::sizeZ; iz++)
            for(int iy=0; iy<FluidBlock::sizeY; iy++)
            for(int ix=0; ix<FluidBlock::sizeX; ix++)
               LHS(ix,iy,iz).s += avgP*h3;
	 }
      }
    }
    else // > 2
    {
       #pragma omp parallel for
       for(size_t i=0; i<vInfo_lhs.size(); ++i)
       {
          ScalarBlock & __restrict__ LHS = *(ScalarBlock*) vInfo_lhs[i].ptrBlock;
          ScalarBlock & __restrict__ Z   = *(ScalarBlock*) vInfo_z  [i].ptrBlock;
          if (vInfo_lhs[i].index[0] == 0 && 
              vInfo_lhs[i].index[1] == 0 && 
              vInfo_lhs[i].index[2] == 0)
          LHS(0,0,0).s = Z(0,0,0).s;
      }
    }
  }
  std::string getName() { return "ComputeLHS"; }
};

class PoissonSolverAMR
{
 protected:
  typedef typename FluidGridMPI::BlockType BlockType;
  SimulationData & sim;
  FluidGridMPI& grid = * sim.grid;
  ComputeLHS findLHS;
  void getZ();
  size_t _dest(const cubism::BlockInfo &info , const int z, const int y, const int x) const
  {
    return BlockType::sizeX * ( BlockType::sizeY * (info.blockID * BlockType::sizeZ  + z) + y) + x;
  }
 public:
  PoissonSolverAMR(SimulationData&s);
  PoissonSolverAMR(const PoissonSolverAMR& c) = delete; 
  void solve();
};

}//namespace cubismup3d
