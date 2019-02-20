//
//  CubismUP_2D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//


#include "PoissonSolver.h"

void PoissonSolver::_cub2fftw() const
{
  const size_t NlocBlocks = local_infos.size();
  #if 0
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
  #endif

  double sumRHS = 0, sumABS = 0;
  #pragma omp parallel for schedule(static) reduction(+ : sumRHS, sumABS)
  for(size_t i=0; i<NlocBlocks; ++i) {
    const size_t offset = _offset(local_infos[i]);
    for(int ix=0; ix<BlockType::sizeX; ix++)
    for(int iy=0; iy<BlockType::sizeY; iy++)
    for(int iz=0; iz<BlockType::sizeZ; iz++) {
      const size_t dest_index = _dest(offset, iz, iy, ix);
      sumABS += std::fabs(data[dest_index]);
      sumRHS +=           data[dest_index];
    }
  }
  double sums[2] = {sumRHS, sumABS};
  MPI_Allreduce(MPI_IN_PLACE, sums, 2, MPI_DOUBLE,MPI_SUM, m_comm);
  sums[1] = std::max(std::numeric_limits<double>::epsilon(), sums[1]);
  const Real correction = sums[0] / sums[1];
  printf("Relative RHS correction:%e / %e\n", sums[0], sums[1]);
  #if 1
    #pragma omp parallel for schedule(static)
    for(size_t i=0; i<NlocBlocks; ++i) {
      const size_t offset = _offset(local_infos[i]);
      for(int ix=0; ix<BlockType::sizeX; ix++)
      for(int iy=0; iy<BlockType::sizeY; iy++)
      for(int iz=0; iz<BlockType::sizeZ; iz++) {
        const size_t dest_index = _dest(offset, iz, iy, ix);
        data[dest_index] -=  std::fabs(data[dest_index]) * correction;
      }
    }

    #ifndef NDEBUG
      double sumRHSpost = 0;
      #pragma omp parallel for schedule(static) reduction(+ : sumRHSpost)
      for(size_t i=0; i<NlocBlocks; ++i) {
        const size_t offset = _offset(local_infos[i]);
        for(int ix=0; ix<BlockType::sizeX; ix++)
        for(int iy=0; iy<BlockType::sizeY; iy++)
        for(int iz=0; iz<BlockType::sizeZ; iz++)
          sumRHSpost += data[_dest(offset, iz, iy, ix)];
      }
      MPI_Allreduce(MPI_IN_PLACE, &sumRHSpost, 1, MPI_DOUBLE,MPI_SUM, m_comm);
      printf("Sum of RHS after correction:%e\n", sumRHSpost);
    #endif
  #endif
}

void PoissonSolver::_fftw2cub() const
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
