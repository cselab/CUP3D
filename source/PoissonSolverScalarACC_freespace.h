//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#pragma once

#include "Definitions.h"
#include "accfft_utils.h"
#include "accfft_common.h"

#include <cuda_runtime_api.h>
#ifndef _FLOAT_PRECISION_
  #include "accfft_gpu.h"
  typedef accfft_plan_gpu myplan;
  typedef Complex myComplex;
  #define accfft_local_size accfft_local_size_dft_r2c_gpu
  #define accfft_plan_dft accfft_plan_dft_3d_r2c_gpu
  #define accfft_clean accfft_cleanup_gpu

#else
  #include "accfft_gpuf.h"
  typedef accfft_plan_gpuf myplan;
  typedef Complexf myComplex;
  #define accfft_local_size accfft_local_size_dft_r2c_gpuf
  #define accfft_plan_dft accfft_plan_dft_3d_r2c_gpuf
  #define accfft_clean accfft_cleanup_gpuf
#endif

using namespace std;

#ifdef _CUDA_COMP_
void _fourier_filter_gpu(myComplex*const __restrict__ data_hat, const int N[3],
  const int isize[3], const int istart[3], const double h);
#endif

template<typename TGrid, typename TStreamer>
class PoissonSolverScalarFFTW_ACC
{
  typedef typename TGrid::BlockType BlockType;
  TGrid& grid;

  const int bs[3] = {BlockType::sizeX, BlockType::sizeY, BlockType::sizeZ};
  const int mybpd[3] = {  grid.getResidentBlocksPerDimension(0),
                          grid.getResidentBlocksPerDimension(1),
                          grid.getResidentBlocksPerDimension(2) };
  const int totbpd[3] = { grid.getBlocksPerDimension(0),
                          grid.getBlocksPerDimension(1),
                          grid.getBlocksPerDimension(2) };
  const size_t myN[3] = { mybpd[0]*bs[0], mybpd[1]*bs[1], mybpd[2]*bs[2] };
  const int totN[3] = { totbpd[0]*bs[0], totbpd[1]*bs[1], totbpd[2]*bs[2] };
  int nprocs, procid, isz[3], osz[3], ist[3], ost[3], alloc_max;
  const vector<BlockInfo> local_infos = grid.getResidentBlocksInfo();
  const size_t N = local_infos.size();
  const int mx = 2 * totN[0] - 1, my = 2 * totN[1] - 1, mz = 2 * totN[2] - 1;
  const int myftNx = (mx+1)/mpisize, mz_pad = (mz/2 +1)*2;
  MPI_Comm c_comm;
  Real * gpuGhat = nullptr;
  Real * cub_rhs = nullptr;
  Real * fft_rhs = nullptr;
  Real * gpu_rhs = nullptr;
  myplan * plan;


public:
  PoissonSolverScalarFFTW_ACC(TGrid& g) : grid(g)
  {
    MPI_Comm_rank(grid.getCartComm(), &procid);
    MPI_Comm_size(grid.getCartComm(), &nprocs);
    int c_dims[2] = { nprocs, 1 };
    accfft_create_comm( grid.getCartComm(), c_dims, &c_comm);
    int M[3] = {mx, my, mz};
    alloc_max = accfft_local_size(M, isz,ist,osz,ost, c_comm);

    printf("[mpi rank %d] isize  %3d %3d %3d osize  %3d %3d %3d\n", procid,
        isz[0],isz[1],isz[2], osz[0],osz[1],osz[2] );
    printf("[mpi rank %d] istart %3d %3d %3d ostart %3d %3d %3d\n", procid,
        ist[0],ist[1],ist[2], ost[0],ost[1],ost[2] );
    assert(isize[0]==(int)myN[0]&&isize[1]==(int)myN[1]&&isize[2]==(int)myN[2]);
    cub_rhs = (Real*)malloc(myN[0] *  myN[1] *  myN[2] * sizeof(Real) );
    fft_rhs = (Real*)malloc(myftNx * totN[1] * totN[2] * sizeof(Real) );
    cudaMalloc((void**) &gpu_rhs, alloc_max);
    cudaMalloc((void**) &gpuGhat, alloc_max / 2);

    plan = accfft_plan_dft(M, gpu_rhs, gpu_rhs, c_comm, ACCFFT_MEASURE);
    initGreen(totN[0],totN[1],totN[2], local_infos[0].h_gridpoint, gpuGhat);
  }

  void solve() {
    _fft2cub(rho);
    dSolveFreespace(*plan, totN[0],totN[1],totN[2], osz[0],osz[1],osz[2],
      gpuGhat, cub_rhs, fft_rhs, gpu_rhs);
  }

  ~PoissonSolverScalarFFTW_ACC() {
    free(cub_rhs);
    free(fft_rhs);
    cudaFree(gpu_rhs);
    cudaFree(gpuGhat);
    accfft_destroy_plan_gpu(plan);
    accfft_clean();
    MPI_Comm_free(&c_comm);
  }

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
    return myIstart[2] + myN[2]*myIstart[1] + myN[2]*myN[1]*myIstart[0];
  }
  inline size_t _dest(const size_t offset,const int z,const int y,const int x) const {
    return offset + z + myN[2] * (y + myN[1] * x);
  }
  inline void _cub2fftw(const size_t offset, const int z, const int y, const int x, const Real rhs) const {
    const size_t dest_index = _dest(offset, z, y, x);
    assert(dest_index >= 0 && dest_index < myN[0] * myN[1] * myN[2]);
    rho[dest_index] = rhs;
  }
  void _fft2cub(Real * out) const {
    #pragma omp parallel for schedule(static)
    for(int i=0; i<N; ++i) {
      const BlockInfo info = local_infos[i];
      BlockType& b = *(BlockType*)info.ptrBlock;
      const size_t offset = _offset(info);
      for(int ix=0; ix<bs[0]; ix++)
      for(int iy=0; iy<bs[1]; iy++)
      for(int iz=0; iz<bs[2]; iz++) {
        const size_t src_index = _dest(offset, iz, iy, ix);
        assert(src_index >= 0 && src_index < myN[0] * myN[1] * myN[2]);
        b(ix,iy,iz).p = out[src_index];
      }
    }
  }
};
