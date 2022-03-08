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
#ifndef CUP_SINGLE_PRECISION
#define MPIREAL MPI_DOUBLE
#else
#define MPIREAL MPI_FLOAT
#endif /* CUP_SINGLE_PRECISION */

namespace cubismup3d {

class ComputeLHS : public Operator
{
  struct KernelLHSPoisson
  {
    const SimulationData & sim;
    KernelLHSPoisson(const SimulationData&s) : sim(s) {}
  
    const StencilInfo stencil{-1,-1,-1,2,2,2,false,{0}};
    void operator()(LabMPIPoisson & lab, const BlockInfo& info, FluidBlockPoisson& o) const
    {
      const double h = info.h; 
      for(int iz=0; iz<FluidBlockPoisson::sizeZ; ++iz)
      for(int iy=0; iy<FluidBlockPoisson::sizeY; ++iy)
      for(int ix=0; ix<FluidBlockPoisson::sizeX; ++ix)
      {
        o(ix,iy,iz).lhs = h*( lab(ix-1,iy,iz).s + lab(ix+1,iy,iz).s + 
                              lab(ix,iy-1,iz).s + lab(ix,iy+1,iz).s +
                              lab(ix,iy,iz-1).s + lab(ix,iy,iz+1).s - 6.0*lab(ix,iy,iz).s);
      }

      BlockCase<FluidBlockPoisson> * tempCase = (BlockCase<FluidBlockPoisson> *)(info.auxiliary);
      FluidBlockPoisson::ElementType * faceXm = nullptr;
      FluidBlockPoisson::ElementType * faceXp = nullptr;
      FluidBlockPoisson::ElementType * faceYm = nullptr;
      FluidBlockPoisson::ElementType * faceYp = nullptr;
      FluidBlockPoisson::ElementType * faceZp = nullptr;
      FluidBlockPoisson::ElementType * faceZm = nullptr;
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
          faceXm[iy + FluidBlock::sizeY * iz].lhs = h*(lab(ix,iy,iz).s - lab(ix-1,iy,iz).s);
        }
      }
      if (faceXp != nullptr)
      {
        int ix = FluidBlock::sizeX-1;
        for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
        for(int iy=0; iy<FluidBlock::sizeY; ++iy)
        {
          faceXp[iy + FluidBlock::sizeY * iz].clear();
          faceXp[iy + FluidBlock::sizeY * iz].lhs = h*(lab(ix,iy,iz).s - lab(ix+1,iy,iz).s);
        }
      }
      if (faceYm != nullptr)
      {
        int iy = 0;
        for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
        for(int ix=0; ix<FluidBlock::sizeX; ++ix)
        {
          faceYm[ix + FluidBlock::sizeX * iz].clear();
          faceYm[ix + FluidBlock::sizeX * iz].lhs = h*(lab(ix,iy,iz).s - lab(ix,iy-1,iz).s);
        }
      }
      if (faceYp != nullptr)
      {
        int iy = FluidBlock::sizeY-1;
        for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
        for(int ix=0; ix<FluidBlock::sizeX; ++ix)
        {
          faceYp[ix + FluidBlock::sizeX * iz].clear();
          faceYp[ix + FluidBlock::sizeX * iz].lhs = h*(lab(ix,iy,iz).s - lab(ix,iy+1,iz).s);
        }
      }
      if (faceZm != nullptr)
      {
        int iz = 0;
        for(int iy=0; iy<FluidBlock::sizeY; ++iy)
        for(int ix=0; ix<FluidBlock::sizeX; ++ix)
        {
          faceZm[ix + FluidBlock::sizeX * iy].clear();
          faceZm[ix + FluidBlock::sizeX * iy].lhs = h*(lab(ix,iy,iz).s - lab(ix,iy,iz-1).s);
        }
      }
      if (faceZp != nullptr)
      {
        int iz = FluidBlock::sizeZ-1;
        for(int iy=0; iy<FluidBlock::sizeY; ++iy)
        for(int ix=0; ix<FluidBlock::sizeX; ++ix)
        {
          faceZp[ix + FluidBlock::sizeX * iy].clear();
          faceZp[ix + FluidBlock::sizeX * iy].lhs = h*(lab(ix,iy,iz).s - lab(ix,iy,iz+1).s);
        }
      }

    }
  };
  public:
  ComputeLHS(SimulationData & s) : Operator(s) { }
  void operator()(const double dt)
  {
    std::vector<cubism::BlockInfo>& vInfoPoisson = gridPoisson->getBlocksInfo();
    const KernelLHSPoisson KPoisson(sim);
    computePoisson<KernelLHSPoisson>(KPoisson,true);

    if (sim.bMeanConstraint == 0) return;

    if (sim.bMeanConstraint <= 2)
    {
       double avgP = 0;
       int index = -1;
       #pragma omp parallel for reduction(+ : avgP)
       for(size_t i=0; i<vInfoPoisson.size(); ++i)
       {
          FluidBlockPoisson & __restrict__ bPoisson  = *(FluidBlockPoisson*) vInfoPoisson[i].ptrBlock;
          const double h3 = vInfoPoisson[i].h*vInfoPoisson[i].h*vInfoPoisson[i].h;
          if (vInfoPoisson[i].index[0] == 0 && 
              vInfoPoisson[i].index[1] == 0 && 
              vInfoPoisson[i].index[2] == 0)
            index = i;
          for(int iz=0; iz<FluidBlock::sizeZ; iz++)
          for(int iy=0; iy<FluidBlock::sizeY; iy++)
          for(int ix=0; ix<FluidBlock::sizeX; ix++)
            avgP += bPoisson(ix,iy,iz).s*h3;
      }
      MPI_Allreduce(MPI_IN_PLACE, &avgP, 1, MPIREAL, MPI_SUM, sim.grid->getCartComm());

      if (sim.bMeanConstraint == 1 && index != -1)
      {
         FluidBlockPoisson & __restrict__ bPoisson  = *(FluidBlockPoisson*) vInfoPoisson[index].ptrBlock;
         bPoisson(0,0,0).lhs = avgP;
      }
      else if (sim.bMeanConstraint == 2)
      {
         #pragma omp parallel for reduction(+ : avgP)
         for(size_t i=0; i<vInfoPoisson.size(); ++i)
	 {
            FluidBlockPoisson & __restrict__ bPoisson  = *(FluidBlockPoisson*) vInfoPoisson[i].ptrBlock;
            const double h3 = vInfoPoisson[i].h*vInfoPoisson[i].h*vInfoPoisson[i].h;
            for(int iz=0; iz<FluidBlock::sizeZ; iz++)
            for(int iy=0; iy<FluidBlock::sizeY; iy++)
            for(int ix=0; ix<FluidBlock::sizeX; ix++)
               bPoisson(ix,iy,iz).lhs += avgP*h3;
	 }
      }
    }
    else // > 2
    {
       #pragma omp parallel for
       for(size_t i=0; i<vInfoPoisson.size(); ++i)
       {
          FluidBlockPoisson & __restrict__ bPoisson  = *(FluidBlockPoisson*) vInfoPoisson[i].ptrBlock;
          if (vInfoPoisson[i].index[0] == 0 && 
              vInfoPoisson[i].index[1] == 0 && 
              vInfoPoisson[i].index[2] == 0)
          bPoisson(0,0,0).lhs = bPoisson(0,0,0).s;
      }
    }
  }
  std::string getName() { return "ComputeLHS"; }
};

class PoissonSolverAMR
{
 protected:
  typedef typename FluidGridMPI::BlockType BlockType;
  typedef typename FluidGridMPIPoisson::BlockType BlockTypePoisson;
  SimulationData & sim;
  FluidGridMPI& grid = * sim.grid;
  FluidGridMPIPoisson& gridPoisson = * sim.gridPoisson;
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
