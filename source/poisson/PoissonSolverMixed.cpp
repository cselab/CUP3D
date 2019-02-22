//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#include "PoissonSolverMixed.h"
#include "PoissonSolver_common.h"


PoissonSolverMixed::PoissonSolverMixed(SimulationData & s) : PoissonSolver(s)
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

  fwd = (void*)_FFTW_(mpi_plan_r2r_3d)(gsize[0], gsize[1], gsize[2], data, data,
    m_comm, XplanF, YplanF, ZplanF, FFTW_MPI_TRANSPOSED_OUT | FFTW_MEASURE);
  bwd = (void*)_FFTW_(mpi_plan_r2r_3d)(gsize[0], gsize[1], gsize[2], data, data,
    m_comm, XplanB, YplanB, ZplanB, FFTW_MPI_TRANSPOSED_IN  | FFTW_MEASURE);

  //std::cout <<    bs[0] << " " <<    bs[1] << " " <<    bs[2] << " ";
  //std::cout <<   myN[0] << " " <<   myN[1] << " " <<   myN[2] << " ";
  //std::cout << gsize[0] << " " << gsize[1] << " " << gsize[2] << " ";
  //std::cout << mybpd[0] << " " << mybpd[1] << " " << mybpd[2] << std::endl;
}

void PoissonSolverMixed::solve()
{
  _cub2fftw();

  _FFTW_(execute)( (fft_plan) fwd);

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

  _FFTW_(execute)( (fft_plan) bwd);

  _fftw2cub();
}

PoissonSolverMixed::~PoissonSolverMixed()
{
  _FFTW_(destroy_plan)((fft_plan) fwd);
  _FFTW_(destroy_plan)((fft_plan) bwd);
  _FFTW_(free)(data);
  _FFTW_(mpi_cleanup)();
}

#undef MPIREAL