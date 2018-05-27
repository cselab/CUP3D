//
//  CubismUP_3D
//
//  Written by Guido Novati ( novatig@ethz.ch ).
//  Copyright (c) 2017 ETHZ. All rights reserved.
//

#pragma once

#include <vector>
#include <cassert>
#include <cmath>
#include <iostream>
#include "Definitions.h"

using namespace std;

#include <BlockInfo.h>
#include <Profiler.h>
#include <pvfmm.hpp>
#include <utils.hpp>


template<typename TGrid, typename TStreamer>
class PoissonSolverScalarPVFMM
{
  typedef typename TGrid::BlockType BlockType;
  TGrid& grid;
  const MPI_Comm comm = grid.getCartComm();
  const int bs[3] = {BlockType::sizeX, BlockType::sizeY, BlockType::sizeZ};
  const vector<BlockInfo> local_infos = grid.getResidentBlocksInfo();
  const size_t N = local_infos.size();
  const size_t mybpd[3] = {
      static_cast<size_t>(grid.getResidentBlocksPerDimension(0)),
      static_cast<size_t>(grid.getResidentBlocksPerDimension(1)),
      static_cast<size_t>(grid.getResidentBlocksPerDimension(2))
  };
  const size_t myN[3]={ mybpd[0]*bs[0], mybpd[1]*bs[1], mybpd[2]*bs[2] };
  const size_t n_local = myN[0]*myN[1]*myN[2], max_pts=300;
  const int mult_order = 8;

  #ifndef _FLOAT_PRECISION_
  const pvfmm::Kernel<double>& kernel_fn = pvfmm::laplace_potn_d;
  #else // _FLOAT_PRECISION_
  const pvfmm::Kernel<float>&  kernel_fn = pvfmm::laplace_potn_f;
  #endif // _FLOAT_PRECISION_

  std::vector<Real>  pts_coord;
  std::vector<Real>  src_value;
  std::vector<Real>  tgt_value;
  std::vector<Real> surf_coord;
  std::vector<Real> surf_value;
  pvfmm::PtFMM_Tree * tree = nullptr;
  // Create memory-manager (not sure what it does, but allocates just 10 Mb)
  pvfmm::mem::MemoryManager mem_mgr(10000000);
  // Load matrices.
  pvfmm::PtFMM matrices(&mem_mgr);

  void setup()
  {
    pts_coord.resize(n_local * 3);
    src_value.resize(n_local, 0);
    tgt_value.resize(n_local, 0);

    #pragma omp parallel for
    for(int i=0; i<N; ++i) {
      const BlockInfo info = local_infos[i];
      const size_t offset = _offset(info);
      for(int iz=0; iz<BlockType::sizeZ; iz++)
      for(int iy=0; iy<BlockType::sizeY; iy++)
      for(int ix=0; ix<BlockType::sizeX; ix++) {
        const size_t pts_index = _dest(offset, iz, iy, ix);
        assert(pts_index>=0 && pts_index<n_local);
        info.pos(&pts_coord[3*pts_index], ix, iy, iz);
      }
    }

    matrices.Initialize(mult_order, comm, &kernel_fn);
  }

  void solve()
  {
    if(tree == nullptr) {
      tree = PtFMM_CreateTree(pts_coord, src_value, surf_coord, surf_value,
                              pts_coord, comm, max_pts, pvfmm::FreeSpace);
      // FMM Setup
      tree->SetupFMM(&matrices);
      PtFMM_Evaluate(tree, tgt_value, n_local);
    }
    else {
      tree->ClearFMMData();
      PtFMM_Evaluate(tree, tgt_value, n_local, &src_value, &surf_value);
    }
    _fftw2cub(data);
  }

  void dispose()
  {
    if(tree not_eq nullptr) delete tree;
  }

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
      info.index[0]*bs[0], info.index[1]*bs[1], info.index[2]*bs[2]
    };
    return myIstart[0] +myN[0]*(myIstart[1] +myN[1]*myIstart[2]);
  }
  inline size_t _dest(const size_t offset,const int z,const int y,const int x) const
  {
    return offset +x +myN[0]*(y +myN[1]*z);
  }
  inline void _cub2fftw(const size_t offset, const int z, const int y, const int x, const Real rhs) const
  {
    const size_t dest_index = _dest(offset, z, y, x);
    assert(dest_index>=0 && dest_index<n_local);
    src_value[dest_index] = rhs;
  }

  void _fftw2cub(const Real * const out) const
  {
    #pragma omp parallel for
    for(int i=0; i<N; ++i) {
      const BlockInfo info = local_infos[i];
      BlockType& b = *(BlockType*)info.ptrBlock;
      const size_t offset = _offset(info);

      for(int iz=0; iz<BlockType::sizeZ; iz++)
      for(int iy=0; iy<BlockType::sizeY; iy++)
      for(int ix=0; ix<BlockType::sizeX; ix++) {
        const size_t dest_index = _dest(offset, iz, iy, ix);
        assert(dest_index>=0 && dest_index<n_local);
        b(ix,iy,iz).p = tgt_value[dest_index];
      }
    }
  }
};
