//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

//#include "PoissonSolverScalarFFTW_ACC.h"
//#include <cuda_runtime_api.h>
#include <array>
#include <vector>
#include <cassert>
#ifndef CUP_SINGLE_PRECISION
  #define MPIREAL MPI_DOUBLE
  typedef double Real;
#else
  #define MPIREAL MPI_FLOAT
  typedef float Real;
#endif
#include "PoissonSolverACC_common.h"


#define CUDA_Check(code) do {  \
    if (code != cudaSuccess) { \
      printf("DONE DEAD func:%s file:%s:%d %s\n", __func__, \
      __FILE__,__LINE__, cudaGetErrorString(code)); \
    } \
} while(0)

__global__
void _fourier_filter_kernel(acc_c*const __restrict__ out,
  const size_t Nx, const size_t Ny, const size_t Nz,
  const size_t nx, const size_t ny, const size_t nz,
  const size_t sx, const size_t sy, const size_t sz,
  const Real wx, const Real wy, const Real wz, const Real fac)
{
  const size_t i = blockDim.x * blockIdx.x + threadIdx.x;
  const size_t j = blockDim.y * blockIdx.y + threadIdx.y;
  const size_t k = blockDim.z * blockIdx.z + threadIdx.z;
  if(i>=nx) return;
  if(j>=ny) return;
  if(k>=nz) return;

  const size_t kx = sx + i, ky = sy + j, kz = sz + k;
  const size_t kkx = kx > Nx/2 ? kx-Nx : kx;
  const size_t kky = ky > Ny/2 ? ky-Ny : ky;
  const size_t kkz = kz > Nz/2 ? kz-Nz : kz;
  const Real rkx = kkx*wx, rky = kky*wy, rkz = kkz*wz;
  const Real kinv = -fac / (rkx*rkx + rky*rky + rkz*rkz);
  //const Real kinv = -scale*(rkx*rkx+rky*rky+rkz*rkz);
  const size_t index = (i*ny + j)*nz + k;
  out[index][0] *= kinv;
  out[index][1] *= kinv;
}


void _fourier_filter_gpu(acc_c*const __restrict__ data_hat,
 const size_t gsize[3],const int osize[3] , const int ostart[3], const double h)
{
  const std::array<Real,3> wfac = {
    Real(2*M_PI/(h*gsize[0])),
    Real(2*M_PI/(h*gsize[1])),
    Real(2*M_PI/(h*gsize[2]))
  };
  const Real scale = 1./( gsize[0]*h * gsize[1]*h * gsize[2]*h );

  int blocksInX = std::ceil(osize[0] / 4.);
  int blocksInY = std::ceil(osize[1] / 4.);
  int blocksInZ = std::ceil(osize[2] / 4.);

  dim3 Dg(blocksInX, blocksInY, blocksInZ);
  dim3 Db(4, 4, 4);
  _fourier_filter_kernel<<<Dg, Db>>>(
    data_hat, gsize[0], gsize[1], gsize[2],
    osize[0], osize[1], osize[2],
    ostart[0], ostart[1], ostart[2],
    wfac[0], wfac[1], wfac[2], scale);

  if(ostart[0]==0 && ostart[1]==0 && ostart[2]==0)
    cudaMemset(data_hat, 0, 2 * sizeof(Real));

  cudaDeviceSynchronize();
}


__global__ void kPos(const int iSzX, const int iSzY, const int iSzZ,
  const int iStX, const int iStY, const int iStZ, const int nGlobX,
  const int nGlobY, const int nGlobZ, const size_t nZpad, Real*const in_out)
{
  const int i = blockDim.x * blockIdx.x + threadIdx.x;
  const int j = blockDim.y * blockIdx.y + threadIdx.y;
  const int k = blockDim.z * blockIdx.z + threadIdx.z;
  if ( (i >= iSzX) || (j >= iSzY) || (k >= iSzZ) ) return;
  const size_t linidx = k + 2*nZpad*(j + iSzY*i);
  const Real I = i + iStX, J = j + iStY, K = k + iStZ;
  in_out[linidx] = K + nGlobZ * (J + nGlobY * I);
}

__global__ void kGreen(const int iSzX, const int iSzY, const int iSzZ,
  const int iStX, const int iStY, const int iStZ,
  const int nGlobX, const int nGlobY, const int nGlobZ, const size_t nZpad,
  const Real h, Real*const in_out)
{
  const int i = blockDim.x * blockIdx.x + threadIdx.x;
  const int j = blockDim.y * blockIdx.y + threadIdx.y;
  const int k = blockDim.z * blockIdx.z + threadIdx.z;
  if ( (i >= iSzX) || (j >= iSzY) || (k >= iSzZ) ) return;
  const size_t linidx = k + 2*nZpad*(j + iSzY*i);
  const int I = i + iStX, J = j + iStY, K = k + iStZ;
  const Real xi = I>=nGlobX? 2*nGlobX-1 - I : I;
  const Real yi = J>=nGlobY? 2*nGlobY-1 - J : J;
  const Real zi = K>=nGlobZ? 2*nGlobZ-1 - K : K;
  const Real r = std::sqrt(xi*xi + yi*yi + zi*zi);
  if(r > 0) in_out[linidx] = - h * h / ( 4 * M_PI * r );
  // G = r_eq^2 / 2 = std::pow(3/8/pi/sqrt(2))^(2/3) * h^2
  else      in_out[linidx] = - Real(0.1924173658) * h * h;
  //else      in_out[linidx] = fac;
}

__global__ void kCopyC2R(const int oSzX,const int oSzY,const int oSzZ,
  const Real norm, const size_t nZpad, const acc_c*const G_hat, Real*const m_kernel)
{
  const int i = threadIdx.x + blockIdx.x * blockDim.x;
  const int j = threadIdx.y + blockIdx.y * blockDim.y;
  const int k = threadIdx.z + blockIdx.z * blockDim.z;
  if ( (i >= oSzX) || (j >= oSzY) || (k >= oSzZ) ) return;
  const size_t linidx = (i*oSzY + j)*nZpad + k;
  m_kernel[linidx] = G_hat[linidx][0] * norm;
}

__global__ void kFreespace(const int oSzX, const int oSzY, const int oSzZ,
  const size_t nZpad, const Real*const G_hat, acc_c*const in_out)
{
  const int i = threadIdx.x + blockIdx.x * blockDim.x;
  const int j = threadIdx.y + blockIdx.y * blockDim.y;
  const int k = threadIdx.z + blockIdx.z * blockDim.z;
  if ( (i >= oSzX) || (j >= oSzY) || (k >= oSzZ) ) return;
  const size_t linidx = (i*oSzY + j)*nZpad + k;
  in_out[linidx][0] *= G_hat[linidx];
  in_out[linidx][1] *= G_hat[linidx];
}

static inline int getSize(MPI_Comm comm) {
  int ret; MPI_Comm_size(comm, &ret); return ret;
}
static inline int getRank(MPI_Comm comm) {
  int ret; MPI_Comm_rank(comm, &ret); return ret;
}

void dSolveFreespace(acc_plan*const P, const int nx,const int ny,const int nz,
  const int locx,const int locy,const int locz, const MPI_Comm comm,
  const int ox,const int oy,const int oz, const Real*const G_hat,
  Real*const cub_rhs, Real*const fft_rhs, Real*const gpu_rhs)
{
  const int mx = 2*nx -1, my = 2*ny -1, mz = 2*nz -1, mz_pad = mz/2 +1;
  const int mpisize = getSize(comm), mpirank = getRank(comm);

  int pos[3], dst[3];
  MPI_Cart_coords(comm, mpirank, 3, pos);
  int szFft[3] = { (mx+1)/mpisize, ny, nz }, start[3]={0,0,0};
  int szCup[3] = { std::min(szFft[0], locx), locy, locz };

  MPI_Datatype submat;
  MPI_Type_create_subarray(3, szFft, szCup, start, MPI_ORDER_C,MPIREAL,&submat);
  MPI_Type_commit(&submat);
  // MPI transfer of data from CUP distribution to 1D-padded FFT distribution
  {
    memset(fft_rhs, 0, szFft[0]*szFft[1]*szFft[2] * sizeof(Real) );
    std::vector<MPI_Request> reqs = std::vector<MPI_Request>(mpisize*2, MPI_REQUEST_NULL);
    const int m_ind =  pos[0]   * locx, m_pos =  mpirank   * szFft[0];
    const int m_nxt = (pos[0]+1)* locx, m_end = (mpirank+1)* szFft[0];
    for(int i=0; i<mpisize; i++)
    {
      MPI_Cart_coords(comm, i, 3, dst); // assert(dst[1]==0 && dst[2]==0);
      const int i_ind =  dst[0]   * locx, i_pos =  i   * szFft[0];
      const int i_nxt = (dst[0]+1)* locx, i_end = (i+1)* szFft[0];
      // test if rank needs to send to i its rhs:
      if( i_pos < m_nxt && m_ind < i_end )
      {
        const int tag = i + mpirank * mpisize;
        const size_t shiftx = std::max(i_pos - m_ind, 0);
        const size_t ptr = szCup[2] * szCup[1] * shiftx;
        MPI_Isend(cub_rhs + ptr, 1, submat, i, tag, comm, &reqs[2*i]);
      }
      // test if rank needs to recv to i's rhs:
      if( m_pos < i_nxt && i_ind < m_end )
      {
        const int tag = mpirank + i * mpisize;
        const size_t shiftx = std::max(i_ind - m_pos, 0);
        const size_t ptr = dst[2]*szCup[2] +nz*(dst[1]*szCup[1] +ny*shiftx);
        MPI_Irecv(fft_rhs + ptr, 1, submat, i, tag, comm, &reqs[2*i + 1]);
      }
    }
    MPI_Waitall(mpisize*2, reqs.data(), MPI_STATUSES_IGNORE);
  }

  // ranks that do not contain only zero-padding, transfer RHS to GPU
  if(mpirank < mpisize/2)
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

  // solve Poisson in padded space
  {
    accfft_exec_r2c(P, gpu_rhs, (acc_c*) gpu_rhs);
    //CUDA_Check(cudaDeviceSynchronize());
    dim3 dB(4, 4, 4);
    dim3 dG(std::ceil(ox/4.), std::ceil(oy/4.), std::ceil(oz/4.));
    kFreespace <<<dG, dB>>> (ox,oy,oz, mz_pad, G_hat, (acc_c*) gpu_rhs);
    //CUDA_Check(cudaDeviceSynchronize());
    accfft_exec_c2r(P, (acc_c*) gpu_rhs, gpu_rhs);
    //CUDA_Check(cudaDeviceSynchronize());
  }

  // ranks that do not contain extended solution, transfer SOL to CPU
  if(mpirank < mpisize/2)
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
  {
    //for(size_t i=0; i<szFft[0]; i++)
    //for(size_t j=0; j<szFft[1]; j++)
    //for(size_t k=0; k<szFft[2]; k++)
    //  fft_rhs[k+szFft[2]*(j+szFft[1]*i)] = k+mz*(j+my*(i+mpirank*szFft[0]));
  }
  // MPI transfer of data from CUP distribution to 1D-padded FFT distribution
  {
    std::vector<MPI_Request> reqs = std::vector<MPI_Request>(mpisize*2, MPI_REQUEST_NULL);
    const int m_ind =  pos[0]   * locx, m_pos =  mpirank   * szFft[0];
    const int m_nxt = (pos[0]+1)* locx, m_end = (mpirank+1)* szFft[0];
    for(int i=0; i<mpisize; i++)
    {
      MPI_Cart_coords(comm, i, 3, dst);
      const int i_ind =  dst[0]   * locx, i_pos =  i   * szFft[0];
      const int i_nxt = (dst[0]+1)* locx, i_end = (i+1)* szFft[0];
      // test if rank needs to send to i its rhs:
      if( i_pos < m_nxt && m_ind < i_end )
      {
        const int tag = i + mpirank * mpisize;
        const size_t shiftx = std::max(i_pos - m_ind, 0);
        const size_t ptr = szCup[2] * szCup[1] * shiftx;
        MPI_Irecv(cub_rhs + ptr, 1, submat, i, tag, comm, &reqs[2*i]);
      }
      // test if rank needs to recv to i's rhs:
      if( m_pos < i_nxt && i_ind < m_end )
      {
        const int tag = mpirank + i * mpisize;
        const size_t shiftx = std::max(i_ind - m_pos, 0);
        const size_t ptr = dst[2]*szCup[2] +nz*(dst[1]*szCup[1] +ny*shiftx);
        MPI_Isend(fft_rhs + ptr, 1, submat, i, tag, comm, &reqs[2*i + 1]);
      }
    }
    MPI_Waitall(mpisize*2, reqs.data(), MPI_STATUSES_IGNORE);
  }
  MPI_Type_free(&submat);
}

void initGreen(const int *isz,const int *osz,const int *ist,const int *ost,
  const int nx,const int ny,const int nz, const Real h, acc_plan* const fwd,
  Real*const m_kernel, Real*const gpu_rhs)
{
  const int mx = 2*nx -1, my = 2*ny -1, mz = 2*nz -1, mz_pad = mz/2 +1;
  {
    dim3 dB(4, 4, 4);
    dim3 dG(std::ceil(isz[0]/4.), std::ceil(isz[1]/4.), std::ceil(isz[2]/4.));
    //cout<<isz[0]<<" "<<isz[1]<<" "<<isz[2]<<" "<< ist[0]<<" "<<ist[1]<<" "<<ist[2]<<" "<<nx<<" "<<ny<<" "<<nz<<" "<<mz_pad<<endl;
    kGreen<<<dG, dB>>> (isz[0],isz[1],isz[2], ist[0],ist[1],ist[2],
      nx, ny, nz, mz_pad, h, gpu_rhs);
    CUDA_Check(cudaDeviceSynchronize());
  }

  accfft_exec_r2c(fwd, gpu_rhs, (acc_c*) gpu_rhs);
  CUDA_Check(cudaDeviceSynchronize());

  {
    const Real norm = 1.0 / (mx*h * my*h * mz*h);
    dim3 dB(4, 4, 4);
    dim3 dG(std::ceil(osz[0]/4.), std::ceil(osz[1]/4.), std::ceil(osz[2]/4.));
    kCopyC2R<<<dG, dB>>> (osz[0],osz[1],osz[2], norm, mz_pad,
      (acc_c*)gpu_rhs, m_kernel);
    CUDA_Check(cudaDeviceSynchronize());
  }
  //{
  //  dim3 dB(4, 4, 4);
  //  dim3 dG(std::ceil(isz[0]/4.), std::ceil(isz[1]/4.), std::ceil(isz[2]/4.));
  //  kPos<<<dG, dB>>> (isz[0],isz[1],isz[2], ist[0],ist[1],ist[2], mx,my,mz, mz_pad, gpu_rhs);
  //}
}
