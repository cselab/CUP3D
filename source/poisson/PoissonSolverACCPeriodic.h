//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#pragma once

#include "PoissonSolver.h"
#include <cuda_runtime_api.h>
#include "PoissonSolverACC_common.h"

void _fourier_filter_gpu(
  acc_c*const __restrict__ data_hat, const size_t gsize[3],
  const int isize[3], const int istart[3], const double h);

class PoissonSolverPeriodic : public PoissonSolver
{
  MPI_Comm c_comm;
  // the local pencil size and the allocation size
  int isize[3], osize[3], istart[3], ostart[3];
  size_t alloc_max;
  Real * rho_gpu;
  Real * phi_gpu;
  acc_c * phi_hat;
  acc_plan * plan;

public:
  PoissonSolverPeriodic(SimulationData & s) : PoissonSolver(s)
  {
    stridez = myN[2];
    stridey = myN[1];

    if (gsize[2]!=myN[2]) {
      printf("Poisson solver assumes grid is distrubuted in x and y directions.\n");
      abort();
    }
    int c_dims[2] = {
      static_cast<int>(gsize[0]/mybpd[0]), static_cast<int>(gsize[1]/mybpd[1])
    };
    assert(gsize[0]%myN[0]==0 && gsize[1]%myN[1]==0);
    accfft_create_comm(grid.getCartComm(), c_dims, &c_comm);
    int _gsize[3] = { static_cast<int>(gsize[0]), static_cast<int>(gsize[1]), static_cast<int>(gsize[2]) };

    alloc_max = accfft_local_size(_gsize, isize, istart, osize, ostart, c_comm);

    //printf("[mpi rank %d] isize  %3d %3d %3d    %3d %3d %3d\n",
    //      m_rank,mybpd[0],mybpd[1],mybpd[2], n[0],n[1],n[2]);
    printf("[mpi rank %d] isize  %3d %3d %3d osize  %3d %3d %3d\n",
      m_rank, isize[0],isize[1],isize[2], osize[0],osize[1],osize[2]
    );
    printf("[mpi rank %d] istart %3d %3d %3d ostart %3d %3d %3d\n",
      m_rank, istart[0],istart[1],istart[2], ostart[0],ostart[1],ostart[2]
    );
    assert(isize[0] == (int) myN[0]);
    assert(isize[1] == (int) myN[1]);
    assert(isize[2] == (int) myN[2]);

    data = (Real*) malloc(isize[0]*isize[1]*isize[2]*sizeof(Real));
    cudaMalloc((void**) &rho_gpu, isize[0]*isize[1]*isize[2]*sizeof(Real));
    //cudaMalloc((void**) &phi_gpu, isize[0]*isize[1]*isize[2]*sizeof(Real));
    cudaMalloc((void**) &phi_hat, alloc_max);

    plan = accfft_plan_dft(_gsize,rho_gpu,(Real*)phi_hat,c_comm,ACCFFT_MEASURE);
  }

  void solve()
  {
    _cub2fftw();

    cudaMemcpy(rho_gpu, data, isize[0]*isize[1]*isize[2]*sizeof(Real),
            cudaMemcpyHostToDevice);

    // Perform forward FFT
    accfft_exec_r2c(plan, rho_gpu, phi_hat);
    // Spectral solve
    _fourier_filter_gpu(phi_hat, gsize, osize, ostart, h);
    // Perform backward FFT
    accfft_exec_c2r(plan, phi_hat, rho_gpu);

    cudaMemcpy(data, rho_gpu, isize[0]*isize[1]*isize[2]*sizeof(Real),
            cudaMemcpyDeviceToHost);

    _fftw2cub();
  }

  ~PoissonSolverPeriodic()
  {
    free(data);
    cudaFree(rho_gpu);
    cudaFree(phi_hat);
    accfft_destroy_plan_gpu(plan);
    accfft_clean();
    MPI_Comm_free(&c_comm);
  }
};

#if 0 // to be maintained tests:
{
  int accfft_left, accfft_right, accfft_bottom, accfft_top;
  //int accfft_front, accfft_back, cubism_front, cubism_back;
  int accfft_rank, accfft_size, cubism_rank;
  MPI_Comm_rank( c_comm, &accfft_rank);
  MPI_Comm_size( c_comm, &accfft_size);
  MPI_Cart_shift(c_comm, 0, 1, &accfft_left,   &accfft_right);
  MPI_Cart_shift(c_comm, 1, 1, &accfft_bottom, &accfft_top);
  //MPI_Cart_shift(c_comm, 2, 1, &accfft_front,  &accfft_back);
  int cubism_left, cubism_right, cubism_bottom, cubism_top;
  MPI_Comm_rank( grid.getCartComm(), &cubism_rank);
  MPI_Cart_shift(grid.getCartComm(), 0, 1, &cubism_left,   &cubism_right);
  MPI_Cart_shift(grid.getCartComm(), 1, 1, &cubism_bottom, &cubism_top);
  //MPI_Cart_shift(grid.getCartComm(), 2, 1, &cubism_front,  &cubism_back);
  //note: accfft comm is not periodic and 2d, cubism is periodic adn 3d, rest must be the same
  if( ( accfft_left  !=MPI_PROC_NULL && accfft_left  !=cubism_left   ) ||
      ( accfft_right !=MPI_PROC_NULL && accfft_right !=cubism_right  ) ||
      ( accfft_bottom!=MPI_PROC_NULL && accfft_bottom!=cubism_bottom ) ||
      ( accfft_top   !=MPI_PROC_NULL && accfft_top   !=cubism_top    ) ||
      ( accfft_rank  !=procid        || accfft_size  !=nprocs        )// ||
    //  ( accfft_front !=MPI_PROC_NULL && accfft_front !=cubism_front  ) ||
    //  ( accfft_back  !=MPI_PROC_NULL && accfft_back  !=cubism_back   )
     ) {
        printf("AccFFT communicator does not match the one from Cubism. Aborting.\n");
        fflush(0);
        MPI_Abort(grid.getCartComm(), MPI_ERR_OTHER);
      }

}
#endif
