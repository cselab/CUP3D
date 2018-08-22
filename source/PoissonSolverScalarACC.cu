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
//}
