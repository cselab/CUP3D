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
#include <array>

#include <cuda_runtime_api.h>
#ifndef CUP_SINGLE_PRECISION
  #include "accfft_gpu.h"
  typedef accfft_plan_gpu myplan;
  #define accfft_local_size accfft_local_size_dft_r2c_gpu
  #define accfft_plan_dft accfft_plan_dft_3d_r2c_gpu
  #define accfft_exec_r2c accfft_execute_r2c_gpu
  typedef accfft_plan_gpu myplan;
  typedef Complex myComplex;
#else
  #include "accfft_gpuf.h"
  typedef accfft_plan_gpuf myplan;
  #define accfft_local_size accfft_local_size_dft_r2c_gpuf
  #define accfft_plan_dft accfft_plan_dft_3d_r2c_gpuf
  #define accfft_exec_r2c accfft_execute_r2c_gpuf
  typedef accfft_plan_gpuf myplan;
  typedef Complexf myComplex;
#endif

typedef Real myCmpl[2];

/******************************************************************************/
void dSolveFreespace(void*const P, const int nx,const int ny,const int nz,
  const int locx,const int locy,const int locz, const MPI_Comm comm,
  const int ox,const int oy,const int oz, const Real*const G_hat,
  Real*const cub_rhs, Real*const fft_rhs, Real*const gpu_rhs);
  void initGreen(const int *isz,const int *osz,const int *ist,const int *ost,
    const int nx,const int ny,const int nz, const Real h, myplan* const fwd,
    Real*const m_kernel, Real*const gpu_rhs);
MPI_Comm my_accfft_create_comm(MPI_Comm C, int c_dims[2]);

void my_accfft_cleanup(void* const plan);

inline void printMemUse(const std::string where) {
  size_t free_byte, total_byte ;
  cudaMemGetInfo( &free_byte, &total_byte ) ;
  double free_db=free_byte, total_db=total_byte, used_db=total_db-free_db;
  printf("%s: used = %f, free = %f MB, total = %f MB\n", where.c_str(),
    used_db/1024/1024, free_db/1024/1024, total_db/1024/1024); fflush(0);
}

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

  const size_t myN[3] = {
    static_cast<size_t>(mybpd[0]*bs[0]),
    static_cast<size_t>(mybpd[1]*bs[1]),
    static_cast<size_t>(mybpd[2]*bs[2])
  };
  const int totN[3] = { totbpd[0]*bs[0], totbpd[1]*bs[1], totbpd[2]*bs[2] };
  const int nprocs = getSize( grid.getCartComm() );
  const int procid = getRank( grid.getCartComm() );

  int isz[3]={0,0,0}, osz[3]={0,0,0}, ist[3]={0,0,0}, ost[3]={0,0,0};
  const std::vector<BlockInfo> local_infos = grid.getResidentBlocksInfo();
  const int mx = 2 * totN[0] - 1, my = 2 * totN[1] - 1, mz = 2 * totN[2] - 1;
  const size_t mz_pad = mz/2 +1, myftNx = (mx+1)/nprocs;

  MPI_Comm c_comm;
  size_t alloc_max;
  Real * cub_rhs;
  Real * fft_rhs;
  Real * gpuGhat;
  Real * gpu_rhs;
  void * plan;

  static int getSize(MPI_Comm comm) {
    int ret; MPI_Comm_size(comm, &ret); return ret;
  }
  static int getRank(MPI_Comm comm) {
    int ret; MPI_Comm_rank(comm, &ret); return ret;
  }
public:
  PoissonSolverScalarFFTW_ACC(TGrid& g) : grid(g)
  {
    int c_dims[2] = { nprocs, 1 };
    c_comm = my_accfft_create_comm( grid.getCartComm(), c_dims);
    int M[3] = {mx, my, mz};
    alloc_max = accfft_local_size(M, isz,ist,osz,ost, c_comm);

    printf("[mpi rank %d] isize  %3d %3d %3d osize  %3d %3d %3d %lu\n", procid,
      isz[0],isz[1],isz[2], osz[0],osz[1],osz[2], alloc_max);
    printf("[mpi rank %d] istart %3d %3d %3d ostart %3d %3d %3d\n", procid,
      ist[0],ist[1],ist[2], ost[0],ost[1],ost[2] );
      fflush(0);

    cudaMalloc((void**) &gpu_rhs, alloc_max);
    cudaMalloc((void**) &gpuGhat, alloc_max/2);
    myplan* fwd = accfft_plan_dft(M, gpu_rhs, gpu_rhs, c_comm, ACCFFT_MEASURE);
    plan = static_cast<void*>(fwd);

    {
      const Real h = grid.getBlocksInfo().front().h_gridpoint;
      initGreen(isz,osz,ist,ost, totN[0],totN[1],totN[2], h, fwd, gpuGhat,gpu_rhs);
    }

    cub_rhs = (Real*) malloc(myN[0]* myN[1]*( myN[2] * sizeof(Real)));
    fft_rhs = (Real*) malloc(myftNx*totN[1]*(totN[2] * sizeof(Real)) );

    if(procid<nprocs/2)
      assert((size_t) isz[0]==myftNx && isz[1]==my && isz[2]==mz);
  }

  void solve()
  {
    const MPI_Comm cart_comm = grid.getCartComm();
    cudaMemset(gpu_rhs, 0, alloc_max);

    #if 0
      Real * cub_test = (Real*) malloc( myN[0]*myN[1]*myN[2] * sizeof(Real) );
      for(size_t i=0;i<myN[0];i++)
        for(size_t j=0;j<myN[1];j++)
          for(size_t k=0;k<myN[2];k++) {
            const auto I = i + procid * myN[0];
            cub_rhs[k+myN[2]*(j+myN[1]*i)]=k+mz*(j+my*I);
          }
      std::copy(cub_rhs, cub_rhs +myN[0]*myN[1]*myN[2], cub_test);
    #endif
    #if 0
      Real * cub_test = (Real*) malloc(myN[0]* myN[1]*(myN[2] * sizeof(Real)));
      std::copy(cub_rhs, cub_rhs +myN[0]*myN[1]*myN[2], cub_test);
      //Real norm = 0;
      //for(size_t i=0; i<myN[0]*myN[1]*myN[2]; i++)
      //  norm+=std::pow(cub_test[i], 2);
      //cout << "norm"<<norm << endl;
    #endif

    dSolveFreespace(plan, totN[0],totN[1],totN[2], myN[0],myN[1],myN[2],
     cart_comm, osz[0],osz[1],osz[2], gpuGhat,cub_rhs,fft_rhs, gpu_rhs);

    #if 0
      Real diff = 0;
      for (size_t i=0; i<myN[0]*myN[1]*myN[2]; i++)
        diff += std::pow(cub_rhs[i]-cub_test[i], 2);
      cout << "diff"<<diff << endl;
      free(cub_test);
    #endif

    _fft2cub();
  }

  ~PoissonSolverScalarFFTW_ACC() {
    free(cub_rhs);
    free(fft_rhs);
    cudaFree(gpu_rhs);
    cudaFree(gpuGhat);
    my_accfft_cleanup(plan);
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
    cub_rhs[dest_index] = rhs;
  }
  void _fft2cub() const {
    #pragma omp parallel for schedule(static)
    for(size_t i=0; i<local_infos.size(); ++i) {
      const BlockInfo info = local_infos[i];
      BlockType& b = *(BlockType*)info.ptrBlock;
      const size_t offset = _offset(info);
      for(int ix=0; ix<bs[0]; ix++)
      for(int iy=0; iy<bs[1]; iy++)
      for(int iz=0; iz<bs[2]; iz++) {
        const size_t src_index = _dest(offset, iz, iy, ix);
        assert(src_index >= 0 && src_index < myN[0] * myN[1] * myN[2]);
        b(ix,iy,iz).p = cub_rhs[src_index];
      }
    }
  }
};
