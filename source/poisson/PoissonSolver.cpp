//
//  CubismUP_2D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//


#include "poisson/PoissonSolver.h"

void PoissonSolver::_cub2fftw() const
{
  assert(stridez>0 && stridey>0 && stridex>0 && data_size>0);
  #if 0
    #pragma omp parallel for schedule(static)
    for(size_t i=0; i<local_infos.size(); ++i) {
      const BlockType& b = *(BlockType*) local_infos[i].ptrBlock;
      const size_t offset = _offset(local_infos[i]);
      for(int iz=0; iz<BlockType::sizeZ; iz++)
      for(int ix=0; ix<BlockType::sizeX; ix++)
      for(int iy=0; iy<BlockType::sizeY; iy++) {
        const size_t dest_index = _dest(offset, iz, iy, ix);
        data[dest_index] = b(ix,iy,iz).p;
      }
    }
  #endif

  Real sumRHS = 0, sumABS = 0;
  #pragma omp parallel for schedule(static) reduction(+ : sumRHS, sumABS)
  for(size_t i=0; i<data_size; ++i) {
    sumABS += std::fabs(data[i]); sumRHS += data[i];
  }
  double sums[2] = {sumRHS, sumABS};
  MPI_Allreduce(MPI_IN_PLACE, sums, 2, MPI_DOUBLE,MPI_SUM, m_comm);
  sums[1] = std::max(std::numeric_limits<double>::epsilon(), sums[1]);
  const Real correction = sums[0] / sums[1];
  if(m_rank == 0)
    printf("Relative RHS correction:%e / %e\n", sums[0], sums[1]);
  #if 1
    #pragma omp parallel for schedule(static)
    for(size_t i=0; i<data_size; ++i) data[i] -= std::fabs(data[i])*correction;

    #ifndef NDEBUG
    {
      double sumRHSpost = 0;
      #pragma omp parallel for schedule(static) reduction(+ : sumRHSpost)
      for(size_t i=0; i<data_size; ++i) sumRHSpost += data[i];
      MPI_Allreduce(MPI_IN_PLACE, &sumRHSpost, 1, MPI_DOUBLE,MPI_SUM, m_comm);
      printf("Sum of RHS after correction:%e\n", sumRHSpost);
    }
    #endif
  #endif
}

void PoissonSolver::_fftw2cub() const
{
  const size_t NlocBlocks = local_infos.size();
  #pragma omp parallel for schedule(static)
  for(size_t i=0; i<NlocBlocks; ++i) {
    BlockType& b = *(BlockType*) local_infos[i].ptrBlock;
    const size_t offset = _offset( local_infos[i] );
    for(int iz=0; iz<BlockType::sizeZ; iz++)
    for(int ix=0; ix<BlockType::sizeX; ix++)
    for(int iy=0; iy<BlockType::sizeY; iy++) {
      const size_t src_index = _dest(offset, iz, iy, ix);
      b(ix,iy,iz).p = data[src_index];
    }
  }
  memset(data, 0, data_size * sizeof(Real));
}
