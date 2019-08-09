//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Hugues de Laroussilhe.
//

#include "SpectralManipACC.h"
#include <cuda_runtime_api.h>
#include "../poisson/PoissonSolverACC_common.h"
#include "accfft_common.h"
#ifndef CUP_SINGLE_PRECISION
  #include "accfft_gpu.h"
  typedef accfft_plan_gpu acc_plan;
#else
  #include "accfft_gpuf.h"
  typedef accfft_plan_gpuf acc_plan;
#endif

#ifndef CUP_SINGLE_PRECISION
#define MPIREAL MPI_DOUBLE
#else
#define MPIREAL MPI_FLOAT
#endif /* CUP_SINGLE_PRECISION */

CubismUP_3D_NAMESPACE_BEGIN
using namespace cubism;

void SpectralManipACC::_compute_largeModesForcing()
{
  Real eps = 0, tke = 0, tkeFiltered = 0;

  // kernel

  MPI_Allreduce(MPI_IN_PLACE, &tkeFiltered, 1, MPIREAL, MPI_SUM, m_comm);
  MPI_Allreduce(MPI_IN_PLACE, &tke, 1, MPIREAL, MPI_SUM, m_comm);
  MPI_Allreduce(MPI_IN_PLACE, &eps, 1, MPIREAL, MPI_SUM, m_comm);

  stats.tke_filtered = tkeFiltered / pow2(normalizeFFT);
  stats.tke = tke / pow2(normalizeFFT);
  stats.eps = eps * 2*(sim.nu) / pow2(normalizeFFT);
}

void SpectralManipACC::_compute_analysis()
{
  Real tke = 0, eps = 0, tauIntegral = 0;
  Real * const E_msr = stats.E_msr;
  memset(E_msr, 0, nBins * sizeof(Real));

  // kernel

  MPI_Allreduce(MPI_IN_PLACE, E_msr, nBins, MPIREAL, MPI_SUM, m_comm);
  MPI_Allreduce(MPI_IN_PLACE, &tke, 1, MPIREAL, MPI_SUM, m_comm);
  MPI_Allreduce(MPI_IN_PLACE, &eps, 1, MPIREAL, MPI_SUM, m_comm);
  MPI_Allreduce(MPI_IN_PLACE, &tauIntegral, 1, MPIREAL, MPI_SUM, m_comm);

  const Real normalization = 1 / pow2(normalizeFFT);
  for (size_t binID = 0; binID < nBins; binID++) E_msr[binID] *= normalization;
  stats.tke = tke * normalization;
  stats.eps = eps * 2 * (sim.nu) * normalization;

  stats.uprime = std::sqrt(2 * stats.tke / 3);
  stats.lambda = std::sqrt(15 * sim.nu / stats.eps) * stats.uprime;
  stats.Re_lambda = stats.uprime * stats.lambda / sim.nu;
  stats.tau_integral = tauIntegral * M_PI/(2*pow3(stats.uprime)) *normalization;
}

void SpectralManipACC::_compute_IC(const std::vector<Real> &K,
                                   const std::vector<Real> &E)
{
  printf("ABORT: SpectralManipACC does not support _compute_IC\n");
  fflush(0); abort(0);
}

SpectralManipACC::SpectralManipACC(SimulationData&s): SpectralManip(s)
{
  if (gsize[2]!=myN[2]) {
    printf("SpectralManipACC assumes grid is distrubuted in x and y.\n");
    abort();
  }
  int c_dims[2] = {
    static_cast<int>(gsize[0]/myN[0]), static_cast<int>(gsize[1]/myN[1])
  };
  assert(gsize[0]%myN[0]==0 && gsize[1]%myN[1]==0);
  accfft_create_comm(grid.getCartComm(), c_dims, &acc_comm);
  int totN[3] = { (int)gsize[0], (int)gsize[1], (int)gsize[2] };

  alloc_max = accfft_local_size(totN, isize, istart, osize, ostart, acc_comm);
  assert(alloc_max == isize[0] * isize[1] * 2*gz_hat * sizeof(Real));

  if(isize[0]!=(int)myN[0] || isize[1]!=(int)myN[1] || isize[2]!=(int)myN[2]) {
    printf("PoissonSolverPeriodic: something wrong in isize\n");
    abort();
  }

  data_size = (size_t) isize[0] * (size_t) isize[1] * (size_t) 2*gz_hat;
  data_u = (Real*) malloc(data_size*sizeof(Real));
  data_v = (Real*) malloc(data_size*sizeof(Real));
  data_w = (Real*) malloc(data_size*sizeof(Real));
  cudaMalloc((void**) & gpu_u, alloc_max);
  cudaMalloc((void**) & gpu_v, alloc_max);
  cudaMalloc((void**) & gpu_w, alloc_max);
  stridez = 1; // fast
  stridey = 2*gz_hat;
  stridex = myN[1] * 2*gz_hat; // slow

  acc_plan* P = accfft_plan_dft(totN, gpu_u, gpu_u, acc_comm, ACCFFT_MEASURE);
  plan = (void*) P;
}

void SpectralManipACC::prepareFwd()
{
}

void SpectralManipACC::prepareBwd()
{
}

void SpectralManipACC::runFwd() const
{
  cudaMemcpy(gpu_u, data_u, alloc_max, cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_v, data_v, alloc_max, cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_w, data_w, alloc_max, cudaMemcpyHostToDevice);
  accfft_exec_r2c((acc_plan*)plan, gpu_u, (acc_c*)gpu_u);
  accfft_exec_r2c((acc_plan*)plan, gpu_v, (acc_c*)gpu_v);
  accfft_exec_r2c((acc_plan*)plan, gpu_w, (acc_c*)gpu_w);
}

void SpectralManipACC::runBwd() const
{
  accfft_exec_c2r((acc_plan*)plan, (acc_c*)gpu_u, gpu_u);
  accfft_exec_c2r((acc_plan*)plan, (acc_c*)gpu_v, gpu_v);
  accfft_exec_c2r((acc_plan*)plan, (acc_c*)gpu_w, gpu_w);
  cudaMemcpy(data_u, gpu_u, alloc_max, cudaMemcpyDeviceToHost);
  cudaMemcpy(data_v, gpu_v, alloc_max, cudaMemcpyDeviceToHost);
  cudaMemcpy(data_w, gpu_w, alloc_max, cudaMemcpyDeviceToHost);
}

SpectralManipACC::~SpectralManipACC()
{
  free(gpu_u);
  free(gpu_v);
  free(gpu_w);
  //cudaFree(rho_gpu);
  cudaFree(gpu_u);
  cudaFree(gpu_v);
  cudaFree(gpu_w);
  accfft_destroy_plan_gpu((acc_plan*)plan);
  accfft_clean();
  MPI_Comm_free(&c_comm);
}

CubismUP_3D_NAMESPACE_END
#undef MPIREAL
