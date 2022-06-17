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
  void getZ();
  size_t _dest(const BlockInfo &info , const int z, const int y, const int x) const
  {
    return BlockType::sizeX * ( BlockType::sizeY * (info.blockID * BlockType::sizeZ  + z) + y) + x;
  }
 public:
  PoissonSolverAMR(SimulationData&s);
  PoissonSolverAMR(const PoissonSolverAMR& c) = delete; 
  void solve();
};

}//namespace cubismup3d
