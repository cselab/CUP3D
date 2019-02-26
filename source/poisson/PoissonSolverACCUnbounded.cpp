//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#include "PoissonSolverACCUnbounded.h"
#include <cuda_runtime_api.h>
#include <cuda_runtime.h>
#include "PoissonSolverACC_common.h"

void dSolveFreespace(const int ox,const int oy,const int oz,const size_t mz_pad,
  const Real*const G_hat, Real*const gpu_rhs);

void initGreen(const int *isz,const int *osz,const int *ist,const int *ost,
    const int nx,const int ny,const int nz, const Real h, acc_plan* const fwd,
    Real*const m_kernel, Real*const gpu_rhs);

PoissonSolverUnbounded::PoissonSolverUnbounded(SimulationData & s) : PoissonSolver(s)
{
  {
    int c_dims[2] = { m_size, 1 };
    accfft_create_comm(m_comm, c_dims, &c_comm);
  }
  {
    int start[3] = {0,0,0};
    MPI_Type_create_subarray(3, szFft,szCup, start,MPI_ORDER_C,MPIREAL,&submat);
    MPI_Type_commit(&submat);
  }
  int M[3] = {mx, my, mz};
  alloc_max = accfft_local_size(M, isize, istart, osize, ostart, c_comm);

  printf("[mpi rank %d] isize  %3d %3d %3d osize  %3d %3d %3d %lu\n", m_rank,
    isize[0],isize[1],isize[2], osize[0],osize[1],osize[2], alloc_max);
  printf("[mpi rank %d] istart %3d %3d %3d ostart %3d %3d %3d\n", m_rank,
    istart[0],istart[1],istart[2], ostart[0],ostart[1],ostart[2] );
    fflush(0);
  //printf("fft sizes %d %d %d cup box %d %d %d\n",szFft[0],szFft[1],szFft[2],szCup[0],szCup[1],szCup[2]);
  cudaMalloc((void**) &gpu_rhs, alloc_max);
  cudaMalloc((void**) &gpuGhat, alloc_max/2);
  acc_plan* P = accfft_plan_dft(M, gpu_rhs, gpu_rhs, c_comm, ACCFFT_MEASURE);
  plan = (void*) P;
  initGreen(isize,osize,istart,ostart, gsize[0],gsize[1],gsize[2], h,P, gpuGhat,gpu_rhs);

  data = (Real*) malloc(myN[0]*  myN[1]*(  myN[2] * sizeof(Real)));
  data_size = (size_t) myN[0] * (size_t) myN[1] * (size_t) myN[2];
  stridez = 1; // fast
  stridey = myN[2];
  stridex = myN[1] * myN[2]; // slow

  fft_rhs = (Real*) malloc(myftNx*gsize[1]*(gsize[2] * sizeof(Real)) );

  if(m_rank<m_size/2)
    assert((size_t) isize[0]==myftNx && isize[1]==my && isize[2]==mz);
}

void PoissonSolverUnbounded::solve()
{
  sim.startProfiler("ACCUNB cub2rhs");
  _cub2fftw();
  cudaMemset(gpu_rhs, 0, alloc_max);
  sim.stopProfiler();

  // MPI transfer of data from CUP distribution to 1D-padded FFT distribution
  sim.startProfiler("ACCUNB rhs2pad");
  cub2padded();
  sim.stopProfiler();

  // ranks that do not contain only zero-padding, transfer RHS to GPU
  sim.startProfiler("ACCUNB cpu2gpu");
  padded2gpu();
  sim.stopProfiler();

  sim.startProfiler("ACCUNB r2c");
  accfft_exec_r2c((acc_plan*) plan, gpu_rhs, (acc_c*) gpu_rhs);
  sim.stopProfiler();

  // solve Poisson in padded Fourier space
  sim.startProfiler("ACCUNB solve");
  dSolveFreespace(osize[0],osize[1],osize[2], mz_pad, gpuGhat, gpu_rhs);
  sim.stopProfiler();

  sim.startProfiler("ACCUNB c2r");
  accfft_exec_c2r((acc_plan*) plan, (acc_c*) gpu_rhs, gpu_rhs);
  sim.stopProfiler();

  // ranks that do not contain extended solution, transfer SOL to CPU
  sim.startProfiler("ACCUNB gpu2cpu");
  gpu2padded();
  sim.stopProfiler();

  sim.startProfiler("ACCUNB pad2rhs");
  padded2cub();
  sim.stopProfiler();

  // MPI transfer of data from 1D-padded FFT distribution to CUP distribution
  sim.startProfiler("ACCUNB rhs2cub");
  _fftw2cub();
  sim.stopProfiler();
}

void PoissonSolverUnbounded::cub2padded() const
{
  int pos[3], dst[3];
  MPI_Cart_coords(m_comm, m_rank, 3, pos);
  memset(fft_rhs, 0, myftNx * gsize[1] * (gsize[2] * sizeof(Real)) );
  std::vector<MPI_Request> reqs = std::vector<MPI_Request>(m_size*2, MPI_REQUEST_NULL);
  const int m_ind =  pos[0]   * myN[0], m_pos =  m_rank   * szFft[0];
  const int m_nxt = (pos[0]+1)* myN[0], m_end = (m_rank+1)* szFft[0];
  for(int i=0; i<m_size; i++)
  {
    MPI_Cart_coords(m_comm, i, 3, dst); // assert(dst[1]==0 && dst[2]==0);
    const int i_ind =  dst[0]   * myN[0], i_pos =  i   * szFft[0];
    const int i_nxt = (dst[0]+1)* myN[0], i_end = (i+1)* szFft[0];
    // test if rank needs to send to i its rhs:
    if( i_pos < m_nxt && m_ind < i_end )
    {
      const int tag = i + m_rank * m_size;
      const size_t shiftx = std::max(i_pos - m_ind, 0);
      //printf("rank %d pos %d %d %d to rank %d pos %d %d %d shift%lu\n", m_rank, pos[0],pos[1],pos[2], i, dst[0],dst[1],dst[2], shiftx); fflush(0);
      const size_t ptr = szCup[2] * szCup[1] * shiftx;
      const size_t num_send = szCup[0] * szCup[1] * szCup[2];
      MPI_Isend(data + ptr, num_send, MPIREAL, i, tag, m_comm, &reqs[2*i]);
    }
    // test if rank needs to recv to i's rhs:
    if( m_pos < i_nxt && i_ind < m_end )
    {
      const int tag = m_rank + i * m_size;
      const size_t shiftx = std::max(i_ind - m_pos, 0);
      //printf("rank %d pos %d %d %d from rank %d pos %d %d %d shift%lu\n", m_rank, pos[0],pos[1],pos[2], i, dst[0],dst[1],dst[2], shiftx); fflush(0);
      const size_t shifty = dst[1]*szCup[1];
      const size_t shiftz = dst[2]*szCup[2];
      const size_t ptr = shiftz + szFft[2]*shifty + szFft[2]*szFft[1]*shiftx;
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
  const int m_ind =  pos[0]   * myN[0], m_pos =  m_rank   * szFft[0];
  const int m_nxt = (pos[0]+1)* myN[0], m_end = (m_rank+1)* szFft[0];
  for(int i=0; i<m_size; i++)
  {
    MPI_Cart_coords(m_comm, i, 3, dst);
    const int i_ind =  dst[0]   * myN[0], i_pos =  i   * szFft[0];
    const int i_nxt = (dst[0]+1)* myN[0], i_end = (i+1)* szFft[0];
    // test if rank needs to send to i its rhs:
    if( i_pos < m_nxt && m_ind < i_end )
    {
      const int tag = i + m_rank * m_size;
      const size_t shiftx = std::max(i_pos - m_ind, 0);
      const size_t ptr = szCup[2] * szCup[1] * shiftx;
      const size_t num_send = szCup[0] * szCup[1] * szCup[2];
      MPI_Irecv(data + ptr, num_send, MPIREAL, i, tag, m_comm, &reqs[2*i]);
    }
    // test if rank needs to recv to i's rhs:
    if( m_pos < i_nxt && i_ind < m_end )
    {
      const int tag = m_rank + i * m_size;
      const size_t shiftx = std::max(i_ind - m_pos, 0);
      const size_t shifty = dst[1]*szCup[1];
      const size_t shiftz = dst[2]*szCup[2];
      const size_t ptr = shiftz + szFft[2]*shifty + szFft[2]*szFft[1]*shiftx;
      MPI_Isend(fft_rhs + ptr, 1, submat, i, tag, m_comm, &reqs[2*i + 1]);
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
  accfft_destroy_plan_gpu((acc_plan*)plan);
  accfft_clean();
  MPI_Comm_free(&c_comm);
  MPI_Type_free(&submat);
}
