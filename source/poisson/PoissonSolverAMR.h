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


#define pVector    tmpU
#define vVector    tmpV
#define rVector    tmpW
#define rhatVector chi
#define xVector    p 
#define sVector    u
#define AxVector   v
#define zVector    w 
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
      const double h = info.h_gridpoint; 
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
    double avgP = 0;
    int index = -1;
    #pragma omp parallel for reduction(+ : avgP)
    for(size_t i=0; i<vInfoPoisson.size(); ++i)
    {
      const bool cornerx = ( vInfoPoisson[i].index[0] == ( (sim.bpdx * (1<<(vInfoPoisson[i].level)) -1)/2 ) ); 
      const bool cornery = ( vInfoPoisson[i].index[1] == ( (sim.bpdy * (1<<(vInfoPoisson[i].level)) -1)/2 ) ); 
      const bool cornerz = ( vInfoPoisson[i].index[2] == ( (sim.bpdz * (1<<(vInfoPoisson[i].level)) -1)/2 ) ); 
      FluidBlockPoisson & __restrict__ bPoisson  = *(FluidBlockPoisson*) vInfoPoisson[i].ptrBlock;
      const double h3 = vInfoPoisson[i].h*vInfoPoisson[i].h*vInfoPoisson[i].h;
      for(int iz=0; iz<FluidBlock::sizeZ; iz++)
      for(int iy=0; iy<FluidBlock::sizeY; iy++)
      for(int ix=0; ix<FluidBlock::sizeX; ix++)
        avgP += bPoisson(ix,iy,iz).s*h3;
      if (cornerx && cornery && cornerz) index = i;
    }
    MPI_Allreduce(MPI_IN_PLACE, &avgP, 1, MPIREAL, MPI_SUM, sim.grid->getWorldComm());
    if (index!=-1)
    {
      FluidBlockPoisson & __restrict__ b  = *(FluidBlockPoisson*) vInfoPoisson[index].ptrBlock;
      b(FluidBlock::sizeX-1,FluidBlock::sizeY-1,FluidBlock::sizeZ-1).lhs = avgP;
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

  const MPI_Comm m_comm = grid.getCartComm();
  const int m_rank = sim.rank, m_size = sim.nprocs;

  Real computeAverage() const;

  ComputeLHS findLHS;
  std::vector<size_t> blocksOffset;
  long long id_min;
  long long id_max;
  size_t iter_min;

 public:
  size_t datasize;

  PoissonSolverAMR(SimulationData&s);
  PoissonSolverAMR(const PoissonSolverAMR& c) = delete; 

  void solve();

  size_t _offset(const cubism::BlockInfo &info) const
  {
    #if 0 //stupid simple and slow approach
      size_t PointsPerBlock = BlockType::sizeX * BlockType::sizeY * BlockType::sizeZ;
      size_t kount = 0;
      std::vector<BlockInfo> & vInfo = grid.getBlocksInfo();
      for (size_t i = 0 ; i < vInfo.size(); i++)
      {
        if (vInfo[i].blockID == info.blockID) break;
        kount ++;
      }
      assert(blocksOffset[info.blockID] == kount * PointsPerBlock);
      return kount * PointsPerBlock;
    #else
      return blocksOffset[info.blockID];
    #endif
  }
  void reset()
  {
    const size_t PointsPerBlock = BlockType::sizeX * BlockType::sizeY * BlockType::sizeZ;
    const size_t Blocks = gridPoisson.getBlocksInfo().size();
    datasize = PointsPerBlock * Blocks;

    std::vector<BlockInfo> & vInfo = gridPoisson.getBlocksInfo();

    id_min = vInfo[0].blockID;
    id_max = vInfo[0].blockID;

    for (size_t i = 1 ; i < vInfo.size(); i++)
    {
      id_min = std::min(vInfo[i].blockID,id_min);
      id_max = std::max(vInfo[i].blockID,id_max);
    }
    blocksOffset.resize(id_max-id_min+1,0);
    for (size_t i = 0 ; i < vInfo.size(); i++)
    {
      blocksOffset[ vInfo[i].blockID-id_min ] = i*PointsPerBlock;
    }
  }
  size_t _dest(const size_t offset,const int z,const int y,const int x) const
  {
    return offset + (BlockType::sizeX*BlockType::sizeY)*z + BlockType::sizeX*y + x;
  }

  void getZ();
};

}//namespace cubismup3d
