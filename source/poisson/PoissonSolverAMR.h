//
//  CubismUP_3D
//  Copyright (c) 2020 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Michalis Chatzimanolakis (michaich@ethz.ch).
//

#pragma once

#include "../SimulationData.h"
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

class PoissonSolverAMR
{
 protected:
  typedef typename FluidGridMPI::BlockType BlockType;
  SimulationData & sim;
  FluidGridMPI& grid = * sim.grid;
  // MPI related
  const MPI_Comm m_comm = grid.getCartComm();
  const int m_rank = sim.rank, m_size = sim.nprocs;

  Real computeAverage() const;
  Real computeRelativeCorrection() const;
 public:
  //Real* data;
  std::vector<Real> data;//initially contains RHS and after solution it contains pressure

  PoissonSolverAMR(SimulationData&s);
  PoissonSolverAMR(const PoissonSolverAMR& c) = delete; 
  ~PoissonSolverAMR();

  void solve();

  size_t _offset(const cubism::BlockInfo &info) const
  {
    assert (info.myrank == m_rank);

    //stupid simple and slow approach, good enough for now
    size_t PointsPerBlock = BlockType::sizeX * BlockType::sizeY * BlockType::sizeZ;
    size_t kount = 0;
    std::vector<BlockInfo> & vInfo = grid.getBlocksInfo();
    for (size_t i = 0 ; i < vInfo.size(); i++)
    {
      if (vInfo[i].blockID == info.blockID) break;
      kount ++;
    }
    return kount * PointsPerBlock;
  }
  void reset()
  {
    const size_t PointsPerBlock = BlockType::sizeX * BlockType::sizeY * BlockType::sizeZ;
    const size_t Blocks = grid.getBlocksInfo().size();
    data.resize(PointsPerBlock * Blocks, 0.0);
  }
  size_t _dest(const size_t offset,const int z,const int y,const int x) const
  {
    return offset + (BlockType::sizeX*BlockType::sizeY)*z + BlockType::sizeX*y + x;
  }
  //void _cub2fftw(const size_t offset, const int z, const int y, const int x, const Real rhs) const
  //{
  //  const size_t dest_index = _dest(offset, z, y, x);
  //  assert(data.size() > dest_index);
  //  data[dest_index] = rhs;
  //}
  Real _fftw2cub(const size_t offset, const int z, const int y, const int x) const
  {
    const size_t dest_index = _dest(offset, z, y, x);
    assert(data.size() > dest_index);
    return data[dest_index];
  }

  void _cub2fftw() const;

  void _fftw2cub() const;


  //will need flux corrections!
  void Get_LHS(bool useX);

  #ifdef PRECOND
  double getA_local(int I1,int I2);
  void FindZ();
  std::vector<std::vector<double>> Ld;
  std::vector <  std::vector <std::vector< std::pair<int,double> > > >L_row;
  std::vector <  std::vector <std::vector< std::pair<int,double> > > >L_col;
  #endif
};

}//namespace cubismup3d