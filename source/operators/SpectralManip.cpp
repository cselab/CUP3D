//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Hugues de Larousssilhe.
//

#include "operators/SpectralManip.h"

#ifndef CUP_SINGLE_PRECISION
#define MPIREAL MPI_DOUBLE
#else
#define MPIREAL MPI_FLOAT
#endif /* CUP_SINGLE_PRECISION */

CubismUP_3D_NAMESPACE_BEGIN
using namespace cubism;

Real energySpectrum::interpE(const Real _k)
{
  int idx_k = -1;
  int size = k.size();
  Real energy = 0.;
  for (int i = 0; i < size; ++i){
    if ( _k < k[i]){
      idx_k = i;
      break;
    }
  }
  if (idx_k > 0)
    energy = E[idx_k] + (E[idx_k-1] - E[idx_k]) / (k[idx_k] - k[idx_k-1]) * (k[idx_k] - _k);
  else if (idx_k==0)
    energy = E[idx_k] - (E[idx_k+1] - E[idx_k]) / (k[idx_k+1] - k[idx_k]) * (k[idx_k+1] - _k);
  else // idx_k==-1
    energy = E[size-1] + (E[size-2] - E[size-1]) / (k[size-1] - k[size-2]) * (k[size-1] - _k);

  energy = (energy < 0) ? 0. : energy;
  return energy;
}

Real energySpectrum::interpSigma2(const Real _k)
{
  int idx_k = -1;
  int size = k.size();
  Real s2 = 0.;
  for (int i = 0; i < size; ++i){
    if ( _k < k[i]){
      idx_k = i;
      break;
    }
  }
  if (idx_k > 0)
    s2 = sigma2[idx_k] + (sigma2[idx_k-1] - sigma2[idx_k]) / (k[idx_k] - k[idx_k-1]) * (k[idx_k] - _k);
  else if (idx_k==0)
    s2 = sigma2[idx_k] - (sigma2[idx_k+1] - sigma2[idx_k]) / (k[idx_k+1] - k[idx_k]) * (k[idx_k] - _k);
  else // idx_k==-1
    s2 = sigma2[size-1] + (sigma2[size-2] - sigma2[size-1]) / (k[size-1] - k[size-2]) * (k[size-1] - _k);

  s2 = (s2 < 0) ? 0. : s2;
  return s2;
}

void energySpectrum::dump2File(const int nBin, const int nGrid, const Real lBox)
{
  const int nBins = std::ceil(std::sqrt(3)*nGrid)/2.0;//*waveFact);
  const Real binSize = M_PI*std::sqrt(3)*nGrid/(nBins*lBox);
  std::stringstream ssR;
  ssR<<"initialSpectrum.dat";
  std::ofstream f;
  f.open (ssR.str());
  f << std::left << std::setw(20) << "k" << std::setw(20) << "Pk" <<std::endl;
  for (int i = 0; i < nBins; ++i){
    const Real k_msr = (i+0.5)*binSize;
    f << std::left << std::setw(20) << k_msr << std::setw(20) << interpE(k_msr) <<std::endl;
  }
}

SpectralManip::SpectralManip(SimulationData & s) : sim(s)
{
  printf("New SpectralManip\n");
  int supported_threads;
  MPI_Query_thread(&supported_threads);
  if (supported_threads<MPI_THREAD_FUNNELED) {
    fprintf(stderr, "SpectralManip ERROR: MPI implementation does not support threads.\n");
    exit(1);
  }

  const int retval = _FFTW_(init_threads)();
  if(retval==0) {
    fprintf(stderr, "SpectralManip: ERROR: Call to fftw_init_threads() returned zero.\n");
    exit(1);
  }
  _FFTW_(mpi_init)();
  const int desired_threads = omp_get_max_threads();
  _FFTW_(plan_with_nthreads)(desired_threads);

  alloc_local = _FFTW_(mpi_local_size_3d_transposed) (
    gsize[0], gsize[1], gsize[2]/2+1, m_comm,
    &local_n0, &local_0_start, &local_n1, &local_1_start);

  data_size = (size_t) myN[0] * (size_t) myN[1] * (size_t) 2*nz_hat;
  stridez = 1; // fast
  stridey = 2*(nz_hat);
  stridex = myN[1] * 2*(nz_hat); // slow

  shifty = static_cast<long>(local_1_start);


  data_u = _FFTW_(alloc_real)(2*alloc_local);
  data_v = _FFTW_(alloc_real)(2*alloc_local);
  data_w = _FFTW_(alloc_real)(2*alloc_local);

  bAllocCs2 = s.bComputeCs2Spectrum;
  if (bAllocCs2)
    data_cs2 = _FFTW_(alloc_real)(2*alloc_local);
}

void SpectralManip::prepareFwd()
{
  if (bAllocFwd) return;

  fwd_u = (void*) _FFTW_(mpi_plan_dft_r2c_3d)(gsize[0], gsize[1], gsize[2],
    data_u, (fft_c*)data_u, m_comm, FFTW_MPI_TRANSPOSED_OUT  | FFTW_MEASURE);
  fwd_v = (void*) _FFTW_(mpi_plan_dft_r2c_3d)(gsize[0], gsize[1], gsize[2],
    data_v, (fft_c*)data_v, m_comm, FFTW_MPI_TRANSPOSED_OUT  | FFTW_MEASURE);
  fwd_w = (void*) _FFTW_(mpi_plan_dft_r2c_3d)(gsize[0], gsize[1], gsize[2],
    data_w, (fft_c*)data_w, m_comm, FFTW_MPI_TRANSPOSED_OUT  | FFTW_MEASURE);

  if (bAllocCs2)
    fwd_cs2 = (void*) _FFTW_(mpi_plan_dft_r2c_3d)(gsize[0], gsize[1], gsize[2],
      data_cs2, (fft_c*)data_cs2, m_comm, FFTW_MPI_TRANSPOSED_OUT  | FFTW_MEASURE);
  bAllocFwd = true;
}

void SpectralManip::prepareBwd()
{
  if (bAllocBwd) return;

  bwd_u = (void*) _FFTW_(mpi_plan_dft_c2r_3d)(gsize[0], gsize[1], gsize[2],
    (fft_c*)data_u, data_u, m_comm, FFTW_MPI_TRANSPOSED_IN  | FFTW_MEASURE);
  bwd_v = (void*) _FFTW_(mpi_plan_dft_c2r_3d)(gsize[0], gsize[1], gsize[2],
    (fft_c*)data_v, data_v, m_comm, FFTW_MPI_TRANSPOSED_IN  | FFTW_MEASURE);
  bwd_w = (void*) _FFTW_(mpi_plan_dft_c2r_3d)(gsize[0], gsize[1], gsize[2],
    (fft_c*)data_w, data_w, m_comm, FFTW_MPI_TRANSPOSED_IN  | FFTW_MEASURE);

  bAllocBwd = true;
}

void SpectralManip::runFwd() const
{
  assert(bAllocFwd);
  // we can use one plan for multiple data:
  //_FFTW_(execute_dft_r2c)( (fft_plan) fwd_u, data_u, (fft_c*)data_u );
  //_FFTW_(execute_dft_r2c)( (fft_plan) fwd_u, data_v, (fft_c*)data_v );
  //_FFTW_(execute_dft_r2c)( (fft_plan) fwd_u, data_w, (fft_c*)data_w );
  _FFTW_(execute)((fft_plan) fwd_u);
  _FFTW_(execute)((fft_plan) fwd_v);
  _FFTW_(execute)((fft_plan) fwd_w);

  if (bAllocCs2)
    _FFTW_(execute)((fft_plan) fwd_cs2);
}

void SpectralManip::runBwd() const
{
  assert(bAllocBwd);
  // we can use one plan for multiple data:
  //_FFTW_(execute_dft_c2r)( (fft_plan) bwd_u, (fft_c*)data_u, data_u );
  //_FFTW_(execute_dft_c2r)( (fft_plan) bwd_u, (fft_c*)data_v, data_v );
  //_FFTW_(execute_dft_c2r)( (fft_plan) bwd_u, (fft_c*)data_w, data_w );
  _FFTW_(execute)((fft_plan) bwd_u);
  _FFTW_(execute)((fft_plan) bwd_v);
  _FFTW_(execute)((fft_plan) bwd_w);
}


SpectralManip::~SpectralManip()
{
  _FFTW_(free)(data_u);
  _FFTW_(free)(data_v);
  _FFTW_(free)(data_w);
  if (bAllocFwd){
  _FFTW_(destroy_plan)((fft_plan) fwd_u);
  _FFTW_(destroy_plan)((fft_plan) fwd_v);
  _FFTW_(destroy_plan)((fft_plan) fwd_w);
  if (bAllocCs2)
    _FFTW_(destroy_plan)((fft_plan) fwd_cs2);
  }
  if (bAllocBwd){
  _FFTW_(destroy_plan)((fft_plan) bwd_u);
  _FFTW_(destroy_plan)((fft_plan) bwd_v);
  _FFTW_(destroy_plan)((fft_plan) bwd_w);
  }

  if (bAllocCs2)
    _FFTW_(free)(data_cs2);

  _FFTW_(mpi_cleanup)();
}
CubismUP_3D_NAMESPACE_END
#undef MPIREAL
