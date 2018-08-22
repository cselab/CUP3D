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
#ifndef _FLOAT_PRECISION_
        typedef double Real;
typedef double Complex[2];
#else
        typedef float Real;
typedef float Complex[2];
#endif

__global__
void _fourier_filter_kernel(Complex*const __restrict__ out,
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

  const Real kinv =(kkx==0&&kky==0&&kkz==0)? 0 : -fac/(rkx*rkx+rky*rky+rkz*rkz);
  //const Real kinv = -scale*(rkx*rkx+rky*rky+rkz*rkz);
  const size_t index = (i*ny + j)*nz + k;
  out[index][0] *= kinv;
  out[index][1] *= kinv;
}


void _fourier_filter_gpu(Complex*const __restrict__ out,
  const std::array<size_t,3> N, const std::array<size_t,3> osize,
  const std::array<size_t,3> ostart, const double h)
{
  const std::array<Real,3> wfac = {
    Real(2*M_PI/(h*N[0])), 
    Real(2*M_PI/(h*N[1])), 
    Real(2*M_PI/(h*N[2]))
  };
  const Real scale = 1./Real(N[0]*N[1]*N[2]);

  int blocksInX = std::ceil(osize[0] / 4.);
  int blocksInY = std::ceil(osize[1] / 4.);
  int blocksInZ = std::ceil(osize[2] / 4.);

  dim3 Dg(blocksInX, blocksInY, blocksInZ);
  dim3 Db(4, 4, 4);
  _fourier_filter_kernel<<<Dg, Db>>>(
    out, N[0], N[1], N[2],
    osize[0], osize[1], osize[2],
    ostart[0], ostart[1], ostart[2],
    wfac[0], wfac[1], wfac[2], scale);

  cudaDeviceSynchronize();
}


__global__ void kGreen(const int iSzX, const int iSzY, const int iSzZ,
  const int iStX, const int iStY, const int iStZ,
  const int nGlobX, const int nGlobY, const int nGlobZ, const int nZpad,
  const const Real fac, const Real h, Real*const in_out) {
  unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;
  unsigned int j = blockDim.y * blockIdx.y + threadIdx.y;
  unsigned int k = blockDim.z * blockIdx.z + threadIdx.z;
  if ( (i >= iSzX) || (j >= iSzY) || (k >= iSzZ) ) return;
  const size_t linidx = k + nZpad*(j + iSzY*i);
  const int I = i + iStX, J = j + iStY, K = k + iStZ;
  const Real xi = I>=nGlobX? 2*nGlobX-1 - I : I;
  const Real yi = J>=nGlobY? 2*nGlobY-1 - J : J;
  const Real zi = K>=nGlobZ? 2*nGlobZ-1 - k : k;
  const Real r = std::sqrt(xi*xi + yi*yi + zi*zi);
  if(r > 0) in_out[linidx] = fac / r;
  // G = r_eq^2 / 2 = std::pow(3/8/pi/sqrt(2))^(2/3) * h^2
  else      in_out[linidx] = 0.1924173658 * h * h;
}

__global__ void kCopyC2R(const int oSzX, const int oSzY, const int oSzZ,
  const Real norm, const Cmpl*const G_hat, Real*const m_kernel) {
  const int i = threadIdx.x + blockIdx.x * blockDim.x;
  const int j = threadIdx.y + blockIdx.y * blockDim.y;
  const int k = threadIdx.z + blockIdx.z * blockDim.z;
  if ( (i >= oSzX) || (j >= oSzY) || (k >= oSzZ) ) return;
  const int linidx = (i*oSzY + j)*oSzZ + k;
  m_kernel[linidx] = G_hat[linidx][0] * norm;
}

__global__ void kFreespace(const int oSzX, const int oSzY, const int oSzZ,
  const Real*const G_hat, Cmpl*const in_out) {
  const int i = threadIdx.x + blockIdx.x * blockDim.x;
  const int j = threadIdx.y + blockIdx.y * blockDim.y;
  const int k = threadIdx.z + blockIdx.z * blockDim.z;
  if ( (i >= oSzX) || (j >= oSzY) || (k >= oSzZ) ) return;
  const int linidx = (i*oSzY + j)*oSzZ + k;
  in_out[linidx][0] *= G_hat[linidx]; in_out[linidx][1] *= G_hat[linidx];
}

void initGreen(const int nx, const int ny, const int nz,
  const Real h, Real*const m_kernel)
{
  const int mx = 2 * nx - 1, my = 2 * ny - 1, mz = 2 * nz - 1;
  const int M[3] = {mx, my, mz}, mz_pad = (mz/2 +1)*2;

  int iSz[3], oSz[3], iSt[3], oSt[3];
  const size_t allocM = accfft_local_size_dft_r2c_gpu(M, iSz,iSt,oSz,oSt, comm);
  cudaMalloc((void**) &tmp, allocM);

  {
    const Real fac = - h * h / ( 4*M_PI );
    dim3 dB(4, 4, 4);
    dim3 dG(std::ceil(iSz[0]/4.), std::ceil(iSz[1]/4.), std::ceil(iSz[2]/4.));
    kGreen<<<dimG, dimB>>> (iSz[0],iSz[1],iSz[2], iSt[0],iSt[1],iSt[2],
      nx, ny, nz, mz_pad, fac, h, tmp);
  }
  {
    myplan fwd;
    fwd = accfft_plan_dft_3d_r2c_gpu(M, tmp, tmp, c_comm, ACCFFT_MEASURE);
    accfft_execute_r2c_gpu(fwd, tmp, (Complex*) tmp);
    accfft_destroy_plan_gpu(fwd);
  }
  {
    const Real norm = 1.0 / (mx * my * mz);
    dim3 dB(4, 4, 4);
    dim3 dG(std::ceil(oSz[0]/4.), std::ceil(oSz[1]/4.), std::ceil(oSz[2]/4.));
    kCopyC2R<<<dG, dB>>> (iSz[0],iSz[1],iSz[2], norm, (Cmpl*)tmp, m_kernel);
  }
  cudaFree(tmp);
}

void dSolveFreespace(const myplan&plan, const int nx,const int ny,const int nz,
  const int ox, const int oy, const int oz, const Real*const G_hat,
  Real*const cub_rhs, Real*const fft_rhs, Real*const gpu_rhs)
{
  const int mx = 2 * nx - 1, my = 2 * ny - 1, mz = 2 * nz - 1;
  const int M[3] = {mx, my, mz}, mz_pad = (mz/2 +1)*2;
  int pos[3], dst[3];
  MPI_Cart_coords(comm, mpirank, 3, pos);
  int szFft[3] = { (mx+1)/mpisize, ny, nz }, start[3]={0,0,0};
  int szCup[3] = { std::min(szFft[0], nx/npx), ny/npy, nz/npz };
  const bool bCast = szFft[0] > nx/npx;
  MPI_Datatype submat;
  MPI_Type_create_subarray(3, szFft, szCup, start, MPI_ORDER_C,MPIREAL,&submat);
  MPI_Type_commit(&submat);
  {
    vector<MPI_Request> reqs = vector<MPI_Request>(mpisize*2, MPI_REQUEST_NULL);
    const int m_ind =  pos[0]   * nx/npx, m_pos =  mpirank   * szFft[0];
    const int m_nxt = (pos[0]+1)* nx/npx, m_end = (mpirank+1)* szFft[0];
    for(int i=0; i<mpisize; i++)
    {
      MPI_Cart_coords(comm, i, 3, dst);
      const int i_ind =  dst[0]   * nx/npx, i_pos =  i   * szFft[0];
      const int i_nxt = (dst[0]+1)* nx/npx, i_end = (i+1)* szFft[0];
      // test if rank needs to send to i its rhs:
      if( i_pos < m_nxt && m_ind < i_end ) {
        const int tag = i + mpirank * mpisize;
        const int shiftx = std::max(i_pos - m_ind, 0);
        const int ptr = szCup[2] * szCup[1] * shiftx;
        MPI_Isend(cub_rhs + ptr, 1, submat, i, tag, comm, &reqs[2*i]);
      }
      // test if rank needs to recv to i's rhs:
      if( m_pos < i_nxt && i_ind < m_end ) {
        const int tag = mpirank + i * mpisize;
        const int shiftx = std::max(i_ind - m_pos, 0);
        const int ptr = dst[2]*szCup[2] +nz*(dst[1]*szCup[1] +ny*shiftx);
        MPI_Irecv(fft_rhs + ptr, 1, submat, i, tag, comm, &reqs[2*i + 1]);
      }
    }
    MPI_Waitall(mpisize*2, reqs.data(), MPI_STATUSES_IGNORE);
  }
  cudaMemset(gpu_rhs, 0, mx*my*mz_pad*sizeof(Real) );
  {
    cudaMemcpy3DParms p = { 0 };
    p.srcPos.x=0; p.srcPos.y=0; p.srcPos.z=0;
    p.dstPos.x=0; p.dstPos.y=0; p.dstPos.z=0;
    p.srcPtr.ptr  = fft_rhs; p.srcPtr.pitch =     nz * sizeof(Real);
    p.srcPtr.xsize = nz;     p.srcPtr.ysize = ny;
    p.dstPtr.ptr  = gpu_rhs; p.dstPtr.pitch = mz_pad * sizeof(Real);
    p.dstPtr.xsize = mz_pad; p.dstPtr.ysize = my;
    p.kind = cudaMemcpyHostToDevice; p.extent.width = nz * sizeof(Real);
    p.extent.height = ny; p.extent.depth = nfftx;
    cudaMemcpy3D(&p)
  }
  {
    accfft_execute_r2c_gpu(plan, gpu_rhs, gpu_rhs);
    dim3 dB(4, 4, 4);
    dim3 dG(std::ceil(iSz[0]/4.), std::ceil(iSz[1]/4.), std::ceil(iSz[2]/4.));
    kFreespace <<<dimG,dimB>>> (ox,oy,oz, G_hat, (Cmpl*)gpu_rhs);
    accfft_execute_c2r_gpu(plan, gpu_rhs, gpu_rhs);
  }
  {
    cudaMemcpy3DParms p = { 0 };
    p.srcPos.x=0; p.srcPos.y=0; p.srcPos.z=0;
    p.dstPos.x=0; p.dstPos.y=0; p.dstPos.z=0;
    p.dstPtr.ptr  = fft_rhs; p.dstPtr.pitch =     nz * sizeof(Real);
    p.dstPtr.xsize = nz;     p.dstPtr.ysize = ny;
    p.srcPtr.ptr  = gpu_rhs; p.srcPtr.pitch = mz_pad * sizeof(Real);
    p.srcPtr.xsize = mz_pad; p.srcPtr.ysize = my;
    p.kind = cudaMemcpyDeviceToHost; p.extent.width = nz * sizeof(Real);
    p.extent.height = ny; p.extent.depth = nfftx;
    cudaMemcpy3D(&p)
  }
  {
    vector<MPI_Request> reqs = vector<MPI_Request>(mpisize*2, MPI_REQUEST_NULL);
    const int m_ind =  pos[0]   * nx/npx, m_pos =  mpirank   * szFft[0];
    const int m_nxt = (pos[0]+1)* nx/npx, m_end = (mpirank+1)* szFft[0];
    for(int i=0; i<mpisize; i++)
    {
      MPI_Cart_coords(comm, i, 3, dst);
      const int i_ind =  dst[0]   * nx/npx, i_pos =  i   * szFft[0];
      const int i_nxt = (dst[0]+1)* nx/npx, i_end = (i+1)* szFft[0];
      // test if rank needs to send to i its rhs:
      if( i_pos < m_nxt && m_ind < i_end ) {
        const int tag = i + mpirank * mpisize;
        const int shiftx = std::max(i_pos - m_ind, 0);
        const int ptr = szCup[2] * szCup[1] * shiftx;
        MPI_Irecv(cub_rhs + ptr, 1, submat, i, tag, comm, &reqs[2*i]);
      }
      // test if rank needs to recv to i's rhs:
      if( m_pos < i_nxt && i_ind < m_end ) {
        const int tag = mpirank + i * mpisize;
        const int shiftx = std::max(i_ind - m_pos, 0);
        const int ptr = dst[2]*szCup[2] +nz*(dst[1]*szCup[1] +ny*shiftx);
        MPI_Isend(fft_rhs + ptr, 1, submat, i, tag, comm, &reqs[2*i + 1]);
      }
    }
    MPI_Waitall(mpisize*2, reqs.data(), MPI_STATUSES_IGNORE);
  }
  MPI_Type_free(&submat);
}
