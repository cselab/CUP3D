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
#define PRECOND

#define pVector  tmpU
#define rVector  tmpV
#define xVector  p //current solution estimate
#define sVector  u
#define AxVector v
#define zVector tmpW //preconditioner


namespace cubismup3d {

class ComputeLHS : public Operator
{
  struct KernelLHS
  {
    const SimulationData & sim;
    KernelLHS(const SimulationData&s) : sim(s) {}
  
    const StencilInfo stencil{-1,-1,-1,2,2,2,false, {FE_TMPU} };
    void operator()(LabMPI & lab, const BlockInfo& info, FluidBlock& o) const
    {
      const double h = info.h_gridpoint; 
      for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
      for(int iy=0; iy<FluidBlock::sizeY; ++iy)
      for(int ix=0; ix<FluidBlock::sizeX; ++ix)
      {
        o(ix,iy,iz).AxVector = h*( lab(ix-1,iy,iz).pVector + lab(ix+1,iy,iz).pVector + 
                                   lab(ix,iy-1,iz).pVector + lab(ix,iy+1,iz).pVector +
                                   lab(ix,iy,iz-1).pVector + lab(ix,iy,iz+1).pVector - 6.0*lab(ix,iy,iz).pVector);
      }

      BlockCase<FluidBlock> * tempCase = (BlockCase<FluidBlock> *)(info.auxiliary);
      FluidBlock::ElementType * faceXm = nullptr;
      FluidBlock::ElementType * faceXp = nullptr;
      FluidBlock::ElementType * faceYm = nullptr;
      FluidBlock::ElementType * faceYp = nullptr;
      FluidBlock::ElementType * faceZp = nullptr;
      FluidBlock::ElementType * faceZm = nullptr;
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
          faceXm[iy + FluidBlock::sizeY * iz].AxVector = h*(lab(ix,iy,iz).pVector - lab(ix-1,iy,iz).pVector);
        }
      }
      if (faceXp != nullptr)
      {
        int ix = FluidBlock::sizeX-1;
        for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
        for(int iy=0; iy<FluidBlock::sizeY; ++iy)
        {
          faceXp[iy + FluidBlock::sizeY * iz].clear();
          faceXp[iy + FluidBlock::sizeY * iz].AxVector = h*(lab(ix,iy,iz).pVector - lab(ix+1,iy,iz).pVector);
        }
      }
      if (faceYm != nullptr)
      {
        int iy = 0;
        for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
        for(int ix=0; ix<FluidBlock::sizeX; ++ix)
        {
          faceYm[ix + FluidBlock::sizeX * iz].clear();
          faceYm[ix + FluidBlock::sizeX * iz].AxVector = h*(lab(ix,iy,iz).pVector - lab(ix,iy-1,iz).pVector);
        }
      }
      if (faceYp != nullptr)
      {
        int iy = FluidBlock::sizeY-1;
        for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
        for(int ix=0; ix<FluidBlock::sizeX; ++ix)
        {
          faceYp[ix + FluidBlock::sizeX * iz].clear();
          faceYp[ix + FluidBlock::sizeX * iz].AxVector = h*(lab(ix,iy,iz).pVector - lab(ix,iy+1,iz).pVector);
        }
      }
      if (faceZm != nullptr)
      {
        int iz = 0;
        for(int iy=0; iy<FluidBlock::sizeY; ++iy)
        for(int ix=0; ix<FluidBlock::sizeX; ++ix)
        {
          faceZm[ix + FluidBlock::sizeX * iy].clear();
          faceZm[ix + FluidBlock::sizeX * iy].AxVector = h*(lab(ix,iy,iz).pVector - lab(ix,iy,iz-1).pVector);
        }
      }
      if (faceZp != nullptr)
      {
        int iz = FluidBlock::sizeZ-1;
        for(int iy=0; iy<FluidBlock::sizeY; ++iy)
        for(int ix=0; ix<FluidBlock::sizeX; ++ix)
        {
          faceZp[ix + FluidBlock::sizeX * iy].clear();
          faceZp[ix + FluidBlock::sizeX * iy].AxVector = h*(lab(ix,iy,iz).pVector - lab(ix,iy,iz+1).pVector);
        }
      }
    }
  };
  public:
  ComputeLHS(SimulationData & s) : Operator(s) { }
  void operator()(const double dt)
  {
    const KernelLHS K(sim);
    compute<KernelLHS>(K,true);
  }
  std::string getName() { return "ComputeLHS"; }
};

class PoissonSolverAMR
{
 protected:
  typedef typename FluidGridMPI::BlockType BlockType;
  SimulationData & sim;
  FluidGridMPI& grid = * sim.grid;

  const MPI_Comm m_comm = grid.getCartComm();
  const int m_rank = sim.rank, m_size = sim.nprocs;

  Real computeAverage() const;
  Real computeRelativeCorrection() const;

  ComputeLHS findLHS;
  std::vector<size_t> blocksOffset;
  long long id_min;
  long long id_max;

 public:
  size_t datasize;

  PoissonSolverAMR(SimulationData&s);
  PoissonSolverAMR(const PoissonSolverAMR& c) = delete; 

  void solve();

  size_t _offset(const cubism::BlockInfo &info) const
  {
    assert (info.myrank == m_rank);
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
    const size_t Blocks = grid.getBlocksInfo().size();
    datasize = PointsPerBlock * Blocks;

    std::vector<BlockInfo> & vInfo = grid.getBlocksInfo();

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

  #ifdef PRECOND
  double getA_local(int I1,int I2);
  void FindZ();
  std::vector<std::vector<double>> Ld;
  std::vector <  std::vector <std::vector< std::pair<int,double> > > >L_row;
  std::vector <  std::vector <std::vector< std::pair<int,double> > > >L_col;
  #endif
};

}//namespace cubismup3d