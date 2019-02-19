//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#pragma once

#include "PoissonSolver.h"

class PoissonSolverMixed : public PoissonSolver
{
  //mycomplex local_rhs, local_work;
  myplan fwd, bwd;
  ptrdiff_t alloc_local=0, local_n0=0, local_0_start=0, local_n1=0, local_1_start=0;

  inline bool DFT_X() const { return sim.BCx_flag == periodic; }
  inline bool DFT_Y() const { return sim.BCy_flag == periodic; }
  inline bool DFT_Z() const { return sim.BCz_flag == periodic; }
 protected:

  template<bool DFTX, bool DFTY, bool DFTZ> void _solve()
  {
    // if BC flag == 1 fourier, else cosine transform
    const Real normX = (DFTX ? 1.0 : 0.5) / ( gsize[0]*h );
    const Real normY = (DFTY ? 1.0 : 0.5) / ( gsize[1]*h );
    const Real normZ = (DFTZ ? 1.0 : 0.5) / ( gsize[2]*h );
    const Real waveFactX = (DFTX ? 2 : 1) * M_PI / ( gsize[0]*h );
    const Real waveFactY = (DFTY ? 2 : 1) * M_PI / ( gsize[1]*h );
    const Real waveFactZ = (DFTZ ? 2 : 1) * M_PI / ( gsize[2]*h );
    const Real norm_factor = normX * normY * normZ;
    Real *const in_out = data;
    const long nKx = static_cast<long>(gsize[0]);
    const long nKy = static_cast<long>(gsize[1]);
    const long nKz = static_cast<long>(gsize[2]);
    const long shifty = static_cast<long>(local_1_start);
    #pragma omp parallel for schedule(static)
    for(long j = 0; j<static_cast<long>(local_n1); ++j)
    for(long i = 0; i<static_cast<long>(gsize[0]); ++i)
    for(long k = 0; k<static_cast<long>(gsize[2]); ++k)
    {
      const size_t linidx = (j*gsize[0] +i)*gsize[2] + k;
      const long J = shifty + j; //memory index plus shift due to decomp
      const long kx = DFTX ? ((i <= nKx/2) ? i : nKx-i) : i;
      const long ky = DFTY ? ((J <= nKy/2) ? J : nKy-J) : J;
      const long kz = DFTZ ? ((k <= nKz/2) ? k : nKz-k) : k;
      const Real rkx = ( kx + (DFTX ? 0 : (Real)0.5 ) ) * waveFactX;
      const Real rky = ( ky + (DFTY ? 0 : (Real)0.5 ) ) * waveFactY;
      const Real rkz = ( kz + (DFTZ ? 0 : (Real)0.5 ) ) * waveFactZ;
      const Real kinv = -1/(rkx*rkx + rky*rky + rkz*rkz);
      in_out[linidx] *= kinv*norm_factor;
    }
    //if (shifty==0 && DFTX && DFTY && DFTZ) in_out[0] = 0;
    if (shifty==0) in_out[0] = 0;
  }

 public:

  PoissonSolverMixed(SimulationData & s) : PoissonSolver(s)
  {
    stridez = myN[2];
    stridey = myN[1];

    int supported_threads;
    MPI_Query_thread(&supported_threads);
    if (supported_threads<MPI_THREAD_FUNNELED) {
    std::cout<<"PoissonSolverMixed ERROR: MPI implementation does not support threads."<<std::endl;
    abort();
    }

    const int retval = _FFTW_(init_threads)();
    if(retval==0) {
      std::cout << "PoissonSolverMixed ERROR: Call to fftw_init_threads() returned zero." << std::endl;
      abort();
    }
    const int desired_threads = omp_get_max_threads();

    _FFTW_(plan_with_nthreads)(desired_threads);
    _FFTW_(mpi_init)();

    alloc_local = _FFTW_(mpi_local_size_3d_transposed) (
      gsize[0], gsize[1], gsize[2], m_comm,
      &local_n0, &local_0_start, &local_n1, &local_1_start);

    auto XplanF = DFT_X() ? FFTW_R2HC : FFTW_REDFT10;
    auto XplanB = DFT_X() ? FFTW_HC2R : FFTW_REDFT01;
    auto YplanF = DFT_Y() ? FFTW_R2HC : FFTW_REDFT10;
    auto YplanB = DFT_Y() ? FFTW_HC2R : FFTW_REDFT01;
    auto ZplanF = DFT_Z() ? FFTW_R2HC : FFTW_REDFT10;
    auto ZplanB = DFT_Z() ? FFTW_HC2R : FFTW_REDFT01;
    data = _FFTW_(alloc_real)(alloc_local);
    fwd = _FFTW_(mpi_plan_r2r_3d)(gsize[0], gsize[1], gsize[2], data, data,
      m_comm, XplanF, YplanF, ZplanF, FFTW_MPI_TRANSPOSED_OUT | FFTW_MEASURE);
    bwd = _FFTW_(mpi_plan_r2r_3d)(gsize[0], gsize[1], gsize[2], data, data,
      m_comm, XplanB, YplanB, ZplanB, FFTW_MPI_TRANSPOSED_IN  | FFTW_MEASURE);

    //std::cout <<    bs[0] << " " <<    bs[1] << " " <<    bs[2] << " ";
    //std::cout <<   myN[0] << " " <<   myN[1] << " " <<   myN[2] << " ";
    //std::cout << gsize[0] << " " << gsize[1] << " " << gsize[2] << " ";
    //std::cout << mybpd[0] << " " << mybpd[1] << " " << mybpd[2] << std::endl;
  }

  void solve()
  {
    _cub2fftw();

    _FFTW_(execute)(fwd);

    if( DFT_X() &&  DFT_Y() &&  DFT_Z()) _solve<1,1,1>();
    else
    if( DFT_X() &&  DFT_Y() && !DFT_Z()) _solve<1,1,0>();
    else
    if( DFT_X() && !DFT_Y() &&  DFT_Z()) _solve<1,0,1>();
    else
    if( DFT_X() && !DFT_Y() && !DFT_Z()) _solve<1,0,0>();
    else
    if(!DFT_X() &&  DFT_Y() &&  DFT_Z()) _solve<0,1,1>();
    else
    if(!DFT_X() &&  DFT_Y() && !DFT_Z()) _solve<0,1,0>();
    else
    if(!DFT_X() && !DFT_Y() &&  DFT_Z()) _solve<0,0,1>();
    else
    if(!DFT_X() && !DFT_Y() && !DFT_Z()) _solve<0,0,0>();
    else {
      printf("Boundary conditions not recognized\n");
      abort();
    }

    _FFTW_(execute)(bwd);

    _fftw2cub();
  }

  ~PoissonSolverMixed()
  {
    _FFTW_(destroy_plan)(fwd);
    _FFTW_(destroy_plan)(bwd);
    _FFTW_(free)(data);
    _FFTW_(mpi_cleanup)();
  }
};
