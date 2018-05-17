//
//  CubismUP_3D
//
//  Written by Guido Novati ( novatig@ethz.ch ).
//  Copyright (c) 2017 ETHZ. All rights reserved.
//

//#include "PoissonSolverScalarFFTW_ACC.h"
//#include <cuda_runtime_api.h>
#ifndef _SP_COMP_
        typedef double Real;
typedef double Complex[2];
#else
        typedef float Real;
typedef float Complex[2];
#endif

__global__
void _fourier_filter_kernel(Complex*const __restrict__ data_hat, const int*const N, const int*const osize, const int*const ostart, const Real*const waveFact, const Real scale)
{
  unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;
  unsigned int j = blockDim.y * blockIdx.y + threadIdx.y;
  unsigned int k = blockDim.z * blockIdx.z + threadIdx.z;
  if(i>=osize[0]) return;
  if(j>=osize[1]) return;
  if(k>=osize[2]) return;

  const int kx = ostart[0]+i;
  const int ky = ostart[1]+j;
  const int kz = ostart[2]+k;
  const int kkx = (kx>N[0]/2) ? kx-N[0] : kx;
  const int kky = (ky>N[1]/2) ? ky-N[1] : ky;
  const int kkz = (kz>N[2]/2) ? kz-N[2] : kz;
  const Real rkx = kkx*waveFact[0];
  const Real rky = kky*waveFact[1];
  const Real rkz = kkz*waveFact[2];

  const Real kinv = (kkx==0&&kky==0&&kkz==0)?0:-scale/(rkx*rkx+rky*rky+rkz*rkz);
  //const Real kinv = -scale*(rkx*rkx+rky*rky+rkz*rkz);
  const long index = (i*osize[1]+j)*osize[2]+k;
  data_hat[index][0] *= kinv;
  data_hat[index][1] *= kinv;
}

void _fourier_filter_gpu(Complex*const __restrict__ data_hat, const int N[3], const int osize[3], const int ostart[3], const double h)
{
  const Real wfac[3] = {2*M_PI/(h*N[0]), 2*M_PI/(h*N[1]), 2*M_PI/(h*N[2])};
  const Real scale = 1./Real(N[0]*N[1]*N[2]);

  int *n_gpu, *osize_gpu, *ostart_gpu;
  Real *fac_gpu;
  cudaMalloc((void**) &n_gpu, 3 * sizeof(int));
  cudaMalloc((void**) &osize_gpu, 3 * sizeof(int));
  cudaMalloc((void**) &ostart_gpu, 3 * sizeof(int));
  cudaMalloc((void**) &fac_gpu, 3 * sizeof(Real));

  cudaMemcpy(n_gpu, N, 3 * sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(osize_gpu, osize, 3 * sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(ostart_gpu, ostart, 3 * sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(fac_gpu, wfac, 3 * sizeof(Real), cudaMemcpyHostToDevice);

  int blocksInX = std::ceil(osize[0] / 4.);
  int blocksInY = std::ceil(osize[1] / 4.);
  int blocksInZ = std::ceil(osize[2] / 4.);
  dim3 Dg(blocksInX, blocksInY, blocksInZ);
  dim3 Db(4, 4, 4);
  _fourier_filter_kernel<<<Dg, Db>>>(data_hat, n_gpu, osize_gpu, ostart_gpu, fac_gpu, scale);

  cudaDeviceSynchronize();
  cudaFree(n_gpu);
  cudaFree(fac_gpu);
  cudaFree(osize_gpu);
  cudaFree(ostart_gpu);
}
//}
