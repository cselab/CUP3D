//
//  CubismUP_3D
//  Copyright (c) 2020 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Michalis Chatzimanolakis (michaich@ethz.ch).
//

#include "PoissonSolverAMR.h"

namespace cubismup3d {

void PoissonSolverAMR::_cub2fftw() const
{
  assert(data.size()>0);
  #ifndef NDEBUG
  {
    computeRelativeCorrection();
  }
  #endif
}

Real PoissonSolverAMR::computeRelativeCorrection() const
{
  Real sumRHS = 0, sumABS = 0;
  const std::vector<cubism::BlockInfo>& vInfo = grid.getBlocksInfo();
  #pragma omp parallel for schedule(static) reduction(+ : sumRHS, sumABS)
  for(size_t i=0; i<vInfo.size(); ++i)
  {
    const double h3 = vInfo[i].h_gridpoint * vInfo[i].h_gridpoint * vInfo[i].h_gridpoint;
    const size_t offset = _offset( vInfo[i] );
    for(int iz=0; iz<BlockType::sizeZ; iz++)
    for(int iy=0; iy<BlockType::sizeY; iy++)
    for(int ix=0; ix<BlockType::sizeX; ix++) {
      const size_t src_index = _dest(offset, iz, iy, ix);
      sumABS += h3 * std::fabs( data[src_index] );
      sumRHS += h3 * data[src_index];
    }
  }
  double sums[2] = {sumRHS, sumABS};
  MPI_Allreduce(MPI_IN_PLACE, sums, 2, MPI_DOUBLE,MPI_SUM, m_comm);
  sums[1] = std::max(std::numeric_limits<double>::epsilon(), sums[1]);
  const Real correction = sums[0] / sums[1];
  if(sim.verbose)
    printf("Relative RHS correction:%e / %e\n", sums[0], sums[1]);
  return correction;
}

Real PoissonSolverAMR::computeAverage() const
{
  Real avgP = 0;
  const std::vector<cubism::BlockInfo>& vInfo = grid.getBlocksInfo();
  #pragma omp parallel for schedule(static) reduction(+ : avgP)
  for(size_t i=0; i<vInfo.size(); ++i)
  {
    const double h3 = vInfo[i].h_gridpoint * vInfo[i].h_gridpoint * vInfo[i].h_gridpoint;
    const size_t offset = _offset( vInfo[i] );
    for(int iz=0; iz<BlockType::sizeZ; iz++)
    for(int iy=0; iy<BlockType::sizeY; iy++)
    for(int ix=0; ix<BlockType::sizeX; ix++) {
      const size_t src_index = _dest(offset, iz, iy, ix);
      avgP += h3 * data[src_index];
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, &avgP, 1, MPIREAL, MPI_SUM, m_comm);
  avgP /= sim.extent[0] * sim.extent[1] * sim.extent[2];
  return avgP;
}

void PoissonSolverAMR::_fftw2cub() const
{
  const std::vector<cubism::BlockInfo>& vInfo = grid.getBlocksInfo();
  #pragma omp parallel for schedule(static)
  for(size_t i=0; i<vInfo.size(); ++i) {
    BlockType& b = *(BlockType*) vInfo[i].ptrBlock;
    const size_t offset = _offset( vInfo[i] );
    for(int iz=0; iz<BlockType::sizeZ; iz++)
    for(int iy=0; iy<BlockType::sizeY; iy++)
    for(int ix=0; ix<BlockType::sizeX; ix++) {
      const size_t src_index = _dest(offset, iz, iy, ix);
      b(ix,iy,iz).p = data[src_index];
    }
  }
}

}//namespace cubismup3d
