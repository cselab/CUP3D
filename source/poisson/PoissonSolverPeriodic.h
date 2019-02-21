//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#pragma once

#include "PoissonSolver.h"
#include "PoissonSolver_common.h"

class PoissonSolverPeriodic : public PoissonSolver
{
  fft_plan fwd, bwd;
  const size_t nz_hat = gsize[2]/2+1;
  const double norm_factor = 1./(gsize[0]*h*gsize[1]*h*gsize[2]*h);
  ptrdiff_t alloc_local=0, local_n0=0, local_0_start=0, local_n1=0, local_1_start=0;

 protected:

  void _solve()
  {
    fft_c *const in_out = (fft_c *) data;
    #if 0
      const Real h2 = h*h;
      const Real factor = h2*norm_factor;

      #pragma omp parallel for
      for(int j = 0; j<local_n1; ++j)
      for(int i = 0; i<gsize[0]; ++i)
      for(int k = 0; k<nz_hat; ++k) {
        const int linidx = (j*gsize[0] +i)*nz_hat + k;
        assert(linidx >=0 && linidx<nx*local_n1*nz_hat);
        assert(linidx < alloc_local);

        const Real denom = 32.*(cos(2.*M_PI*i/nx) + cos(2.*M_PI*(local_1_start+j)/ny) + cos(2.*M_PI*k/nz)) -
                  2.*(cos(4.*M_PI*i/nx) + cos(4.*M_PI*(local_1_start+j)/ny) + cos(4.*M_PI*k/nz)) - 90.;

        const Real inv_denom = (denom==0)? 0.:1./denom;
        const Real fatfactor = 12. * inv_denom * factor;

        in_out[linidx][0] *= fatfactor;
        in_out[linidx][1] *= fatfactor;
      }
    #else
      const Real waveFactX = 2.0*M_PI/(gsize[0]*h);
      const Real waveFactY = 2.0*M_PI/(gsize[1]*h);
      const Real waveFactZ = 2.0*M_PI/(gsize[2]*h);
      const long nKx = static_cast<long>(gsize[0]);
      const long nKy = static_cast<long>(gsize[1]);
      const long nKz = static_cast<long>(gsize[2]);
      const long shifty = static_cast<long>(local_1_start);
      #pragma omp parallel for
      for(long j = 0; j<static_cast<long>(local_n1); ++j)
      for(long i = 0; i<static_cast<long>(gsize[0]); ++i)
      for(long k = 0; k<static_cast<long>(nz_hat);   ++k) {
        const size_t linidx = (j*gsize[0] +i)*nz_hat + k;
        const long kx = (i <= nKx/2) ? i : -(nKx-i);
        const long l = shifty + j; //memory index plus shift due to decomp
        const long ky = (l <= nKy/2) ? l : -(nKy-l);
        const long kz = (k <= nKz/2) ? k : -(nKz-k);

        const Real rkx = kx*waveFactX, rky = ky*waveFactY, rkz = kz*waveFactZ;
        const Real kinv =  -1/(rkx*rkx+rky*rky+rkz*rkz);
        in_out[linidx][0] *= kinv*norm_factor;
        in_out[linidx][1] *= kinv*norm_factor;
      }
    #endif

    //this is sparta!
    if (local_1_start == 0) in_out[0][0] = in_out[0][1] = 0;
  }

 public:

  PoissonSolverPeriodic(SimulationData & s) : PoissonSolver(s)
  {
    stridez = 2*nz_hat;
    stridey = myN[1];

    int supported_threads;
    MPI_Query_thread(&supported_threads);
    if (supported_threads<MPI_THREAD_FUNNELED) {
      std::cout<<"PoissonSolverPeriodic ERROR: MPI implementation does not support threads."<<std::endl;
      abort();
    }

    const int retval = _FFTW_(init_threads)();
    if(retval==0) {
      std::cout << "PoissonSolverPeriodic: ERROR: Call to fftw_init_threads() returned zero." << std::endl;
      abort();
    }
    const int desired_threads = omp_get_max_threads();
    _FFTW_(plan_with_nthreads)(desired_threads);
    _FFTW_(mpi_init)();

    alloc_local = _FFTW_(mpi_local_size_3d_transposed) (
      gsize[0], gsize[1], gsize[2]/2+1, m_comm,
      &local_n0, &local_0_start, &local_n1, &local_1_start);

    data = _FFTW_(alloc_real)(2*alloc_local);
    fwd = _FFTW_(mpi_plan_dft_r2c_3d)(gsize[0], gsize[1], gsize[2],
      data, (fft_c *)data, m_comm, FFTW_MPI_TRANSPOSED_OUT | FFTW_MEASURE);
    bwd = _FFTW_(mpi_plan_dft_c2r_3d)(gsize[0], gsize[1], gsize[2],
      (fft_c *)data, data, m_comm, FFTW_MPI_TRANSPOSED_IN  | FFTW_MEASURE);

    //std::cout <<    bs[0] << " " <<    bs[1] << " " <<    bs[2] << " ";
    //std::cout <<   myN[0] << " " <<   myN[1] << " " <<   myN[2] << " ";
    //std::cout << gsize[0] << " " << gsize[1] << " " << gsize[2] << " ";
    //std::cout << mybpd[0] << " " << mybpd[1] << " " << mybpd[2] << std::endl;
  }

  void solve()
  {
    _cub2fftw();

    _FFTW_(execute)(fwd);

    _solve();

    _FFTW_(execute)(bwd);

    _fftw2cub();
  }

  ~PoissonSolverPeriodic()
  {
    _FFTW_(destroy_plan)(fwd);
    _FFTW_(destroy_plan)(bwd);
    _FFTW_(free)(data);
    _FFTW_(mpi_cleanup)();
  }
};
