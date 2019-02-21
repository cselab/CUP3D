//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#pragma once

#include "PoissonSolverACCUnbounded.h"
#include <cuda_runtime_api.h>
#include "PoissonSolverACC_common.h"

inline void printMemUse(const std::string where)
{
  size_t free_byte, total_byte ;
  cudaMemGetInfo( &free_byte, &total_byte ) ;
  double free_db=free_byte, total_db=total_byte, used_db=total_db-free_db;
  printf("%s: used = %f, free = %f MB, total = %f MB\n", where.c_str(),
    used_db/1024/1024, free_db/1024/1024, total_db/1024/1024); fflush(0);
}

void dSolveFreespace(const int ox,const int oy,const int oz,const size_t mz_pad,
  const Real*const G_hat, acc_c*const gpu_rhs);

void initGreen(const int *isz,const int *osz,const int *ist,const int *ost,
    const int nx,const int ny,const int nz, const Real h, acc_plan* const fwd,
    Real*const m_kernel, Real*const gpu_rhs);

PoissonSolverUnbounded::PoissonSolverUnbounded(SimulationData & s) : PoissonSolver(s)
{
  stridez = myN[2];
  stridey = myN[1];
  int c_dims[2] = { m_size, 1 };
  accfft_create_comm(m_comm, c_dims, &c_comm);
  MPI_Type_create_subarray(3, szFft, szCup, start, MPI_ORDER_C,MPIREAL,&submat);
  MPI_Type_commit(&submat);
  int M[3] = {mx, my, mz};
  alloc_max = accfft_local_size(M, isize, istart, osize, ostart, c_comm);

  printf("[mpi rank %d] isize  %3d %3d %3d osize  %3d %3d %3d %lu\n", m_rank,
    isize[0],isize[1],isize[2], osize[0],osize[1],osize[2], alloc_max);
  printf("[mpi rank %d] istart %3d %3d %3d ostart %3d %3d %3d\n", m_rank,
    istart[0],istart[1],istart[2], ostart[0],ostart[1],ostart[2] );
    fflush(0);

  cudaMalloc((void**) &gpu_rhs, alloc_max);
  cudaMalloc((void**) &gpuGhat, alloc_max/2);
  acc_plan* P = accfft_plan_dft(M, gpu_rhs, gpu_rhs, c_comm, ACCFFT_MEASURE);
  plan = (void*) P;
  initGreen(isize,osize,istart,ostart, gsize[0],gsize[1],gsize[2], h, plan, gpuGhat, gpu_rhs);

  data = (Real*) malloc(myN[0]*  myN[1]*(  myN[2] * sizeof(Real)));
  fft_rhs = (Real*) malloc(myftNx*gsize[1]*(gsize[2] * sizeof(Real)) );

  if(m_rank<m_size/2)
    assert((size_t) isize[0]==myftNx && isize[1]==my && isize[2]==mz);
}

void PoissonSolverUnbounded::solve()
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
  // MPI transfer of data from CUP distribution to 1D-padded FFT distribution
  cub2padded();
  // ranks that do not contain only zero-padding, transfer RHS to GPU
  padded2gpu();
  accfft_exec_r2c((acc_plan*) plan, gpu_rhs, (acc_c*) gpu_rhs);
  // solve Poisson in padded Fourier space
  dSolveFreespace(osize[0],osize[1],osize[2], mz_pad, gpuGhat, gpu_rhs);
  accfft_exec_c2r((acc_plan*) plan, (acc_c*) gpu_rhs, gpu_rhs);
  // ranks that do not contain extended solution, transfer SOL to CPU
  gpu2padded();
  // MPI transfer of data from 1D-padded FFT distribution to CUP distribution
  padded2cub();
  #if 0
    Real diff = 0;
    for (size_t i=0; i<myN[0]*myN[1]*myN[2]; i++)
      diff += std::pow(data[i]-cub_test[i], 2);
    cout << "diff"<<diff << endl;
    free(cub_test);
  #endif

  _fftw2cub();
}

void PoissonSolverUnbounded::cub2padded() const
{
  int pos[3], dst[3];
  MPI_Cart_coords(m_comm, m_rank, 3, pos);

  memset(fft_rhs, 0, szFft[0]*szFft[1]*szFft[2] * sizeof(Real) );
  std::vector<MPI_Request> reqs = std::vector<MPI_Request>(m_size*2, MPI_REQUEST_NULL);
  const int m_ind =  pos[0]   * locx, m_pos =  m_rank   * szFft[0];
  const int m_nxt = (pos[0]+1)* locx, m_end = (m_rank+1)* szFft[0];
  for(int i=0; i<m_size; i++)
  {
    MPI_Cart_coords(m_comm, i, 3, dst); // assert(dst[1]==0 && dst[2]==0);
    const int i_ind =  dst[0]   * locx, i_pos =  i   * szFft[0];
    const int i_nxt = (dst[0]+1)* locx, i_end = (i+1)* szFft[0];
    // test if rank needs to send to i its rhs:
    if( i_pos < m_nxt && m_ind < i_end )
    {
      const int tag = i + m_rank * m_size;
      const size_t shiftx = std::max(i_pos - m_ind, 0);
      const size_t ptr = szCup[2] * szCup[1] * shiftx;
      MPI_Isend(data + ptr, 1, submat, i, tag, m_comm, &reqs[2*i]);
    }
    // test if rank needs to recv to i's rhs:
    if( m_pos < i_nxt && i_ind < m_end )
    {
      const int tag = m_rank + i * m_size;
      const size_t shiftx = std::max(i_ind - m_pos, 0);
      const size_t ptr = dst[2]*szCup[2] +nz*(dst[1]*szCup[1] +ny*shiftx);
      MPI_Irecv(fft_rhs + ptr, 1, submat, i, tag, m_comm, &reqs[2*i + 1]);
    }
  }
  MPI_Waitall(m_size*2, reqs.data(), MPI_STATUSES_IGNORE);
}

void PoissonSolverUnbounded::padded2cub() const
{
  int pos[3], dst[3];
  MPI_Cart_coords(m_comm, m_rank, 3, pos);

  std::vector<MPI_Request> reqs = std::vector<MPI_Request>(m_size*2, MPI_REQUEST_NULL);
  const int m_ind =  pos[0]   * locx, m_pos =  m_rank   * szFft[0];
  const int m_nxt = (pos[0]+1)* locx, m_end = (m_rank+1)* szFft[0];
  for(int i=0; i<m_size; i++)
  {
    MPI_Cart_coords(comm, i, 3, dst);
    const int i_ind =  dst[0]   * locx, i_pos =  i   * szFft[0];
    const int i_nxt = (dst[0]+1)* locx, i_end = (i+1)* szFft[0];
    // test if rank needs to send to i its rhs:
    if( i_pos < m_nxt && m_ind < i_end )
    {
      const int tag = i + m_rank * m_size;
      const size_t shiftx = std::max(i_pos - m_ind, 0);
      const size_t ptr = szCup[2] * szCup[1] * shiftx;
      MPI_Irecv(data + ptr, 1, submat, i, tag, comm, &reqs[2*i]);
    }
    // test if rank needs to recv to i's rhs:
    if( m_pos < i_nxt && i_ind < m_end )
    {
      const int tag = m_rank + i * m_size;
      const size_t shiftx = std::max(i_ind - m_pos, 0);
      const size_t ptr = dst[2]*szCup[2] +nz*(dst[1]*szCup[1] +ny*shiftx);
      MPI_Isend(fft_rhs + ptr, 1, submat, i, tag, comm, &reqs[2*i + 1]);
    }
  }
  MPI_Waitall(m_size*2, reqs.data(), MPI_STATUSES_IGNORE);
}

void PoissonSolverUnbounded::padded2gpu() const
{
  if(m_rank < m_size/2)
  {
  #if 1
    cudaMemcpy3DParms p = {};
    p.srcPos.x=0; p.srcPos.y=0; p.srcPos.z=0; p.dstPos.x=0; p.dstPos.y=0; p.dstPos.z=0;
    p.dstPtr = make_cudaPitchedPtr(gpu_rhs, 2*mz_pad*sizeof(Real), 2*mz_pad, my);
    p.srcPtr = make_cudaPitchedPtr(fft_rhs, szFft[2]*sizeof(Real), szFft[2], szFft[1]);
    p.extent = make_cudaExtent(szFft[2]*sizeof(Real), szFft[1], szFft[0]);
    p.kind = cudaMemcpyHostToDevice;
    CUDA_Check(cudaMemcpy3D(&p));
  #else
    for(int i=0; i<szFft[0]; i++) {
      CUDA_Check(cudaMemcpy2D(
        gpu_rhs + 2*mz_pad*my*i, 2*mz_pad*sizeof(Real),
        fft_rhs + szFft[2]*szFft[1]*i, szFft[2]*sizeof(Real),
        szFft[2]*sizeof(Real), szFft[1], // sizes
        cudaMemcpyHostToDevice) );
    }
  #endif
    //CUDA_Check(cudaDeviceSynchronize());
  }
}

void PoissonSolverUnbounded::gpu2padded() const
{
  if(m_rank < m_size/2)
  {
  #if 1
    cudaMemcpy3DParms p = {};
    p.srcPos.x=0; p.srcPos.y=0; p.srcPos.z=0; p.dstPos.x=0; p.dstPos.y=0; p.dstPos.z=0;
    p.srcPtr = make_cudaPitchedPtr(gpu_rhs, 2*mz_pad*sizeof(Real), 2*mz_pad, my);
    p.dstPtr = make_cudaPitchedPtr(fft_rhs, szFft[2]*sizeof(Real), szFft[2], szFft[1]);
    p.extent = make_cudaExtent(szFft[2]*sizeof(Real), szFft[1], szFft[0]);
    p.kind = cudaMemcpyDeviceToHost;
    CUDA_Check(cudaMemcpy3D(&p));
  #else
    for(int i=0; i<szFft[0]; i++) {
      CUDA_Check(cudaMemcpy2D(
        fft_rhs + szFft[2]*szFft[1]*i, szFft[2]*sizeof(Real),
        gpu_rhs + 2*mz_pad*my*i, 2*mz_pad*sizeof(Real),
        szFft[2]*sizeof(Real), szFft[1], // sizes
        cudaMemcpyDeviceToHost) );
    }
  #endif
    //CUDA_Check(cudaDeviceSynchronize());
  }
}
PoissonSolverUnbounded::~PoissonSolverUnbounded()
{
  free(data);
  free(fft_rhs);
  cudaFree(gpu_rhs);
  cudaFree(gpuGhat);
  accfft_destroy_plan_gpu(plan);
  accfft_clean();
  MPI_Comm_free(&c_comm);
  MPI_Type_free(&submat);
}
