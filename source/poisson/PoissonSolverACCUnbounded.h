//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#pragma once

#include "Definitions.h"
#include <array>

#include <cuda_runtime_api.h>


/******************************************************************************/
void dSolveFreespace(acc_plan*const P, const int nx,const int ny,const int nz,
  const int locx,const int locy,const int locz, const MPI_Comm comm,
  const int ox,const int oy,const int oz, const Real*const G_hat,
  Real*const cub_rhs, Real*const fft_rhs, Real*const gpu_rhs);

void initGreen(const int *isz,const int *osz,const int *ist,const int *ost,
    const int nx,const int ny,const int nz, const Real h, acc_plan* const fwd,
    Real*const m_kernel, Real*const gpu_rhs);

inline void printMemUse(const std::string where)
{
  size_t free_byte, total_byte ;
  cudaMemGetInfo( &free_byte, &total_byte ) ;
  double free_db=free_byte, total_db=total_byte, used_db=total_db-free_db;
  printf("%s: used = %f, free = %f MB, total = %f MB\n", where.c_str(),
    used_db/1024/1024, free_db/1024/1024, total_db/1024/1024); fflush(0);
}

class PoissonSolverUnbounded : public PoissonSolver
{
  // the local pencil size and the allocation size
  int isize[3], osize[3], istart[3], ostart[3];
  const int mx = 2*gsize[0]-1, my = 2*gsize[1]-1, mz = 2*gsize[2]-1;
  const size_t mz_pad = mz/2 +1, myftNx = (mx+1)/m_size;

  MPI_Comm c_comm;
  size_t alloc_max;
  Real * fft_rhs;
  Real * gpuGhat;
  Real * gpu_rhs;
  acc_plan * plan;

public:
  PoissonSolverUnbounded(SimulationData & s) : PoissonSolver(s)
  {
    stridez = myN[2];
    stridey = myN[1];
    int c_dims[2] = { m_size, 1 };
    accfft_create_comm(m_comm, c_dims, &c_comm);
    int M[3] = {mx, my, mz};
    alloc_max = accfft_local_size(M, isize, istart, osize, ostart, c_comm);

    printf("[mpi rank %d] isize  %3d %3d %3d osize  %3d %3d %3d %lu\n", m_rank,
      isize[0],isize[1],isize[2], osize[0],osize[1],osize[2], alloc_max);
    printf("[mpi rank %d] istart %3d %3d %3d ostart %3d %3d %3d\n", m_rank,
      istart[0],istart[1],istart[2], ostart[0],ostart[1],ostart[2] );
      fflush(0);

    cudaMalloc((void**) &gpu_rhs, alloc_max);
    cudaMalloc((void**) &gpuGhat, alloc_max/2);
    plan = accfft_plan_dft(M, gpu_rhs, gpu_rhs, c_comm, ACCFFT_MEASURE);

    initGreen(isize,osize,istart,ostart, gsize[0],gsize[1],gsize[2], h, plan, gpuGhat, gpu_rhs);

    data = (Real*) malloc(myN[0]*  myN[1]*(  myN[2] * sizeof(Real)));
    fft_rhs = (Real*) malloc(myftNx*gsize[1]*(gsize[2] * sizeof(Real)) );

    if(m_rank<m_size/2)
      assert((size_t) isize[0]==myftNx && isize[1]==my && isize[2]==mz);
  }

  void solve()
  {
    _cub2fftw();

    cudaMemset(gpu_rhs, 0, alloc_max);

    #if 0
      Real * cub_test = (Real*) malloc( myN[0]*myN[1]*myN[2] * sizeof(Real) );
      for(size_t i=0;i<myN[0];i++)
        for(size_t j=0;j<myN[1];j++)
          for(size_t k=0;k<myN[2];k++) {
            const auto I = i + procid * myN[0];
            data[k+myN[2]*(j+myN[1]*i)]=k+mz*(j+my*I);
          }
      std::copy(data, data +myN[0]*myN[1]*myN[2], cub_test);
    #endif
    #if 0
      Real * cub_test = (Real*) malloc(myN[0]* myN[1]*(myN[2] * sizeof(Real)));
      std::copy(data, data +myN[0]*myN[1]*myN[2], cub_test);
      //Real norm = 0;
      //for(size_t i=0; i<myN[0]*myN[1]*myN[2]; i++)
      //  norm+=std::pow(cub_test[i], 2);
      //cout << "norm"<<norm << endl;
    #endif

    dSolveFreespace(plan, gsize[0],gsize[1],gsize[2], myN[0],myN[1],myN[2],
     m_comm, osize[0],osize[1],osize[2], gpuGhat, data, fft_rhs, gpu_rhs);

    #if 0
      Real diff = 0;
      for (size_t i=0; i<myN[0]*myN[1]*myN[2]; i++)
        diff += std::pow(data[i]-cub_test[i], 2);
      cout << "diff"<<diff << endl;
      free(cub_test);
    #endif

    _fftw2cub();
  }

  ~PoissonSolverUnbounded()
  {
    free(data);
    free(fft_rhs);
    cudaFree(gpu_rhs);
    cudaFree(gpuGhat);
    accfft_destroy_plan_gpu(plan);
    accfft_clean();
    MPI_Comm_free(&c_comm);
  }
};
