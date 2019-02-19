//
//  CubismUP_2D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#pragma once

#include <vector>
#include <cassert>
#include <cstring>

#include <fftw3-mpi.h>

#include "../SimulationData.h"

#include "Cubism/BlockInfo.h"

#ifndef CUP_SINGLE_PRECISION
#define _FFTW_(s) fftw_##s
typedef fftw_complex mycomplex;
typedef fftw_plan myplan;
#define MPIREAL MPI_DOUBLE
#else
#define _FFTW_(s) fftwf_##s
typedef fftwf_complex mycomplex;
typedef fftwf_plan myplan;
#define MPIREAL MPI_FLOAT
#endif /* CUP_SINGLE_PRECISION */

class PoissonSolver
{
 protected:
  typedef typename FluidGridMPI::BlockType BlockType;
  const SimulationData & sim;
  FluidGridMPI& grid = * sim.grid;
  // MPI related
  const MPI_Comm m_comm;
  const int m_rank, m_size;

  static constexpr int bs[3] = {BlockType::sizeX, BlockType::sizeY, BlockType::sizeZ};
  const std::vector<BlockInfo> local_infos = grid.getResidentBlocksInfo();

  const size_t mybpd[3] = {
      static_cast<size_t>(grid.getResidentBlocksPerDimension(0)),
      static_cast<size_t>(grid.getResidentBlocksPerDimension(1)),
      static_cast<size_t>(grid.getResidentBlocksPerDimension(2))
  };
  const size_t gsize[3] = {
      static_cast<size_t>(grid.getBlocksPerDimension(0)*bs[0]),
      static_cast<size_t>(grid.getBlocksPerDimension(1)*bs[1]),
      static_cast<size_t>(grid.getBlocksPerDimension(2)*bs[2])
  };
  const size_t myN[3]={ mybpd[0]*bs[0], mybpd[1]*bs[1], mybpd[2]*bs[2] };
  const double h = grid.getBlocksInfo().front().h_gridpoint;

  size_t stridez = 0;
  size_t stridey = 0;
  Real* data;

 public:
  PoissonSolver(SimulationData&s) : sim(s), m_comm(grid.getCartComm()),
  m_rank(s.rank), m_size(s.nprocs)
  {
    if (StreamerDiv::channels != 1) {
      std::cout << "PoissonSolverScalar_MPI(): Error: StreamerDiv::channels is "
                << StreamerDiv::channels << " (should be 1)." << std::endl;
      abort();
    }
  }
  PoissonSolver(const PoissonSolver& c) = delete;
  virtual ~PoissonSolver() {}

  virtual void solve() = 0;

  inline size_t _offset(const int blockID) const
  {
    const BlockInfo &info = local_infos[blockID];
    return _offset(info);
  }
  inline size_t _offset_ext(const BlockInfo &info) const
  {
    assert(local_infos[info.blockID].blockID == info.blockID);
    return _offset(local_infos[info.blockID]);
    //for(const auto & local_info : local_infos)
    //  if(local_info.blockID == info.blockID)
    //    return _offset(local_info);
    //printf("PSolver cannot find obstacle block\n");
    //abort();
  }
  inline size_t _offset(const BlockInfo &info) const
  {
    const int myIstart[3] = {
      info.index[0]*bs[0],
      info.index[1]*bs[1],
      info.index[2]*bs[2]
    };
    return myIstart[2] +stridez*(myIstart[1]+ stridey*myIstart[0]);
  }
  inline size_t _dest(const size_t offset,const int z,const int y,const int x) const
  {
    return offset + z + stridez *(y + stridey * x);
  }

  inline void _cub2fftw(const size_t offset, const int z, const int y, const int x, const Real rhs) const
  {
    const size_t dest_index = _dest(offset, z, y, x);
    data[dest_index] = rhs;
  }
  inline Real _fftw2cub(const size_t offset, const int z, const int y, const int x) const
  {
    const size_t dest_index = _dest(offset, z, y, x);
    return data[dest_index];
  }

  void _cub2fftw() const
  {
    const size_t NlocBlocks = local_infos.size();
    #pragma omp parallel for
    for(size_t i=0; i<NlocBlocks; ++i) {
      const BlockInfo info = local_infos[i];
      BlockType& b = *(BlockType*)info.ptrBlock;
      const size_t offset = _offset(info);

      for(int ix=0; ix<BlockType::sizeX; ix++)
      for(int iy=0; iy<BlockType::sizeY; iy++)
      for(int iz=0; iz<BlockType::sizeZ; iz++) {
        const size_t dest_index = _dest(offset, iz, iy, ix);
        data[dest_index] = b(ix,iy,iz).p;
      }
    }
  }
  void _fftw2cub() const
  {
    const size_t NlocBlocks = local_infos.size();
    #pragma omp parallel for
    for(size_t i=0; i<NlocBlocks; ++i) {
      const BlockInfo info = local_infos[i];
      BlockType& b = *(BlockType*)info.ptrBlock;
      const size_t offset = _offset(info);

      for(int ix=0; ix<BlockType::sizeX; ix++)
      for(int iy=0; iy<BlockType::sizeY; iy++)
      for(int iz=0; iz<BlockType::sizeZ; iz++) {
        const size_t src_index = _dest(offset, iz, iy, ix);
        b(ix,iy,iz).p = data[src_index];
      }
    }
  }
  //  assert(src_index>=0 && src_index<gsize[0]*gsize[1]*gsize[2]);
  //  assert(dest_index>=0 && dest_index<gsize[0]*gsize[1]*nz_hat*2);
  // assert(dest_index < m_local_N0*m_NN1*2*m_Nzhat);
};
