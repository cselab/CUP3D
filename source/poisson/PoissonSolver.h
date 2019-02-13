//
//  CubismUP_2D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//


#pragma once

#include "../Operator.h"

class PoissonSolver
{
 protected:
  SimulationData& sim;

  // total number of DOFs in y and x:
  const size_t totNz = sim.vel->getBlocksPerDimension(2) * VectorBlock::sizeZ;
  const size_t totNy = sim.vel->getBlocksPerDimension(1) * VectorBlock::sizeY;
  const size_t totNx = sim.vel->getBlocksPerDimension(0) * VectorBlock::sizeX;
  const size_t stride;

  // memory buffer for mem transfers to/from solver:
  Real * buffer = nullptr; // rhs in cub2rhs, sol in sol2cub
  Real * const presMom = new Real[totNy * totNx];
  void cub2rhs();
  void sol2cub();

 public:
  int iter = 0;
  PoissonSolver(SimulationData& s, long stride);

  virtual void solve() = 0;

  virtual ~PoissonSolver() { }

  inline size_t _offset(const int blockID) const {
    const BlockInfo &info = local_infos[blockID];
    return _offset(info);
  }
  inline size_t _offset_ext(const BlockInfo &info) const {
    assert(local_infos[info.blockID].blockID == info.blockID);
    return _offset(local_infos[info.blockID]);
    //for(const auto & local_info : local_infos)
    //  if(local_info.blockID == info.blockID)
    //    return _offset(local_info);
    //printf("PSolver cannot find obstacle block\n");
    //abort();
  }
  inline size_t _offset(const BlockInfo &info) const {
    const int myIstart[3] = {
      info.index[0]*bs[0],
      info.index[1]*bs[1],
      info.index[2]*bs[2]
    };
    return myIstart[2] +2*nz_hat*(myIstart[1]+ myN[1]*myIstart[0]);
  }
  inline size_t _dest(const size_t offset,const int z,const int y,const int x) const
  {
    return offset + z + 2*nz_hat*(y + myN[1] * x);
  }
  inline void _cub2fftw(const size_t offset, const int z, const int y, const int x, const Real rhs) const
  {
    const size_t dest_index = _dest(offset, z, y, x);
    assert(dest_index>=0 && dest_index<gsize[0]*gsize[1]*nz_hat*2);
    data[dest_index] = rhs;
  }

  void _cub2fftw(Real * const out) const
  {
    #pragma omp parallel for
    for(size_t i=0; i<local_infos.size(); ++i) {
      const BlockInfo info = local_infos[i];
      BlockType& b = *(BlockType*)info.ptrBlock;
      const size_t offset = _offset(info);

      for(int ix=0; ix<BlockType::sizeX; ix++)
      for(int iy=0; iy<BlockType::sizeY; iy++)
      for(int iz=0; iz<BlockType::sizeZ; iz++) {
        const size_t dest_index = _dest(offset, iz, iy, ix);
        assert(dest_index>=0 && dest_index<gsize[0]*gsize[1]*nz_hat*2);
        out[dest_index] = b(ix,iy,iz).p;
        //TStreamer::operate(b.data[iz][iy][ix], &out[dest_index]);
      }
    }
  }

  inline size_t _offset(const int blockID) const
  {
      const BlockInfo &info = m_local_infos[blockID];
      return _offset(info);
  }

  inline size_t _offset_ext(const BlockInfo &info) const
  {
      assert(m_local_infos[info.blockID].blockID == info.blockID);
      return _offset(m_local_infos[info.blockID]);
  }

  inline size_t _offset(const BlockInfo &info) const
  {
      // info must be a local BlockInfo! (obtained from
      // grid->getResidentBlocksInfo())
      const int myIstart[3] = {
          info.index[0] * BlockType::sizeX,
          info.index[1] * BlockType::sizeY,
          info.index[2] * BlockType::sizeZ
      };
      return m_fft.get_data_linear_offset(myIstart[0], myIstart[1], myIstart[2]);
  }

  inline size_t _dest(const size_t offset,const int z,const int y,const int x) const
  {
      return offset + m_fft.get_data_linear_offset(x, y, z);
  }

  inline void _cub2fftw(const size_t offset, const int z, const int y, const int x, const Real rhs)
  {
      const size_t dest_index = _dest(offset, z, y, x);
      m_fft.put_data(dest_index, rhs);
  }

  void _cub2fftw()
  {
    #pragma omp parallel for
      for(size_t i=0; i<m_local_infos.size(); ++i)
      {
          const BlockInfo info = m_local_infos[i];
          BlockType& b = *(BlockType*)info.ptrBlock;
          const size_t offset = _offset(info);

          for(int ix=0; ix<BlockType::sizeX; ix++)
              for(int iy=0; iy<BlockType::sizeY; iy++)
                  for(int iz=0; iz<BlockType::sizeZ; iz++)
                  {
                      const size_t dest_index = _dest(offset, iz, iy, ix);
                      m_fft.put_data(dest_index, b(ix,iy,iz).p);
                  }
      }
  }
};
