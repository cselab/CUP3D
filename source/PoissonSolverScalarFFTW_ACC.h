//
//  CubismUP_3D
//
//  Written by Guido Novati ( novatig@ethz.ch ).
//  Copyright (c) 2017 ETHZ. All rights reserved.
//

#pragma once
//#include <mpi.h>
//#include <vector>
//#include <cassert>
//#include <cmath>
//#include <iostream>

#include "Definitions.h"
#include "accfft_utils.h"
#include "accfft_common.h"

#ifdef _CUDA_COMP_
#include <cuda_runtime_api.h>
  #ifndef _SP_COMP_
    #include "accfft_gpu.h"
    typedef accfft_plan_gpu myplan;
    typedef Complex myComplex;
  #else
    #include "accfft_gpuf.h"
    typedef accfft_plan_gpuf myplan;
    typedef Complexf myComplex;
  #endif
#else
  #ifndef _SP_COMP_
    #include "accfft.h"
    typedef accfft_plan myplan;
    typedef Complex myComplex;
  #else
    #include "accfftf.h"
    typedef accfft_planf myplan;
    typedef Complexf myComplex;
  #endif
#endif

using namespace std;

#ifdef _CUDA_COMP_
void _fourier_filter_gpu(myComplex*const __restrict__ data_hat, const int N[3],
  const int isize[3], const int istart[3], const double h);
#endif

template<typename TGrid, typename TStreamer>
class PoissonSolverScalarFFTW_ACC
{
  typedef typename TGrid::BlockType BlockType;
  TGrid& grid;

  const int bs[3], mybpd[3], totbpd[3];
  const size_t myN[3]={ mybpd[0]*bs[0], mybpd[1]*bs[1], mybpd[2]*bs[2] };

  int nprocs, procid, isize[3],osize[3],istart[3],ostart[3], alloc_max;
  const vector<BlockInfo> local_infos = grid.getResidentBlocksInfo();
  const size_t N = local_infos.size();

  MPI_Comm c_comm;
  Real * rho;
  #ifdef _CUDA_COMP_
  Real * rho_gpu;
  Real * phi_gpu;
  #endif
  myComplex * phi_hat;
  myplan * plan;

  void _cub2fft(Real * out) const
  {
    #pragma omp parallel for
    for(int i=0; i<N; ++i) {
      const BlockInfo info = local_infos[i];
      BlockType& b = *(BlockType*)info.ptrBlock;
      const size_t offset = _offset(info);

      for(int ix=0; ix<bs[0]; ix++)
      for(int iy=0; iy<bs[1]; iy++)
      for(int iz=0; iz<bs[2]; iz++) {
        const size_t dest_index = _dest(offset, iz, iy, ix);
        assert(dest_index >= 0 && dest_index < myN[0] * myN[1] * myN[2]);
        TStreamer::operate(b.data[iz][iy][ix], &out[dest_index]);
      }
    }
  }

  void _fft2cub(Real * out) const
  {
    #pragma omp parallel for
    for(int i=0; i<N; ++i) {
      const BlockInfo info = local_infos[i];
      BlockType& b = *(BlockType*)info.ptrBlock;
      const size_t offset = _offset(info);

      for(int ix=0; ix<bs[0]; ix++)
      for(int iy=0; iy<bs[1]; iy++)
      for(int iz=0; iz<bs[2]; iz++) {
        const size_t src_index = _dest(offset, iz, iy, ix);
        assert(src_index >= 0 && src_index < myN[0] * myN[1] * myN[2]);
        //b(ix,iy,iz).p = factor * out[src_index];
        b(ix,iy,iz).p = out[src_index];
      }
    }
  }

  #ifndef _CUDA_COMP_
  void _fourier_filter(myComplex*const __restrict__ data_hat, const Real h)
  {
    const int NX = totbpd[0]*bs[0];
    const int NY = totbpd[1]*bs[1];
    const int NZ = totbpd[2]*bs[2];
    const Real waveFactX = 2.0*M_PI/(h*NX);
    const Real waveFactY = 2.0*M_PI/(h*NY);
    const Real waveFactZ = 2.0*M_PI/(h*NZ);
    const Real norm_factor = 1./Real(NX*NY*NZ);

    #pragma omp parallel for collapse(3)
    for (int i=0; i < osize[0]; i++)
    for (int j=0; j < osize[1]; j++)
    for (int k=0; k < osize[2]; k++) {
      const int kx = ostart[0]+i;
      const int ky = ostart[1]+j;
      const int kz = ostart[2]+k;
      const int kkx = (kx>NX/2) ? kx-NX : kx;
      const int kky = (ky>NY/2) ? ky-NY : ky;
      const int kkz = (kz>NZ/2) ? kz-NZ : kz;
      const Real rkx = kkx*waveFactX;
      const Real rky = kky*waveFactY;
      const Real rkz = kkz*waveFactZ;
      const Real kinv = (kkx==0 && kky==0 && kkz==0) ? 0.0
              : -norm_factor/(rkx*rkx+rky*rky+rkz*rkz);
      const int index = (i*osize[1]+j)*osize[2]+k;
      data_hat[index][0] *= kinv;
      data_hat[index][1] *= kinv;
    }
  }
  #endif

public:
  PoissonSolverScalarFFTW_ACC(const int desired_threads, TGrid& g)
  : grid(g), bs{BlockType::sizeX, BlockType::sizeY, BlockType::sizeZ},
    mybpd{g.getResidentBlocksPerDimension(0), g.getResidentBlocksPerDimension(1), g.getResidentBlocksPerDimension(2)},
    totbpd{g.getBlocksPerDimension(0), g.getBlocksPerDimension(1), g.getBlocksPerDimension(2)}
  {
    if (totbpd[2]!=mybpd[2]) {
      printf("Poisson solver assumes grid is distrubuted in x and y directions.\n");
      abort();
    }

    MPI_Comm_rank(grid.getCartComm(), &procid);
    MPI_Comm_size(grid.getCartComm(), &nprocs);
    int n[3] = {totbpd[0]*bs[0], totbpd[1]*bs[1], totbpd[2]*bs[2]};
    int loc[3] = {mybpd[0]*bs[0], mybpd[1]*bs[1], mybpd[2]*bs[2]};
    int c_dims[2] = { totbpd[0]/mybpd[0], totbpd[1]/mybpd[1] };
    assert(totbpd[0]%mybpd[0]==0 && totbpd[1]%mybpd[1]==0);
    accfft_create_comm(grid.getCartComm(),c_dims,&c_comm);
    {
      int accfft_left, accfft_right, accfft_bottom, accfft_top, accfft_front, accfft_back, accfft_rank, accfft_size;
      MPI_Comm_rank( c_comm, &accfft_rank);
      MPI_Comm_size( c_comm, &accfft_size);
      MPI_Cart_shift(c_comm, 0, 1, &accfft_left,   &accfft_right);
      MPI_Cart_shift(c_comm, 1, 1, &accfft_bottom, &accfft_top);
      //MPI_Cart_shift(c_comm, 2, 1, &accfft_front,  &accfft_back);
      int cubism_left, cubism_right, cubism_bottom, cubism_top, cubism_front, cubism_back, cubism_rank;
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

    // Get the local pencil size and the allocation size
    #ifdef _CUDA_COMP_
      #ifndef _SP_COMP_
      alloc_max = accfft_local_size_dft_r2c_gpu(n,isize,istart,osize,ostart,c_comm);
      #else
      alloc_max = accfft_local_size_dft_r2c_gpuf(n,isize,istart, osize, ostart,c_comm);
      #endif
    #else
      #ifndef _SP_COMP_
      alloc_max = accfft_local_size_dft_r2c(n,isize,istart,osize,ostart,c_comm);
      #else
      alloc_max = accfft_local_size_dft_r2cf(n,isize,istart,osize,ostart,c_comm);
      #endif
    #endif
    //printf("[mpi rank %d] isize  %3d %3d %3d    %3d %3d %3d\n",
    //      procid,mybpd[0],mybpd[1],mybpd[2], n[0],n[1],n[2]);
    printf("[mpi rank %d] isize  %3d %3d %3d osize  %3d %3d %3d\n", procid,
        isize[0],isize[1],isize[2],
        osize[0],osize[1],osize[2]
    );
    printf("[mpi rank %d] istart %3d %3d %3d ostart %3d %3d %3d\n", procid,
        istart[0],istart[1],istart[2],
        ostart[0],ostart[1],ostart[2]
    );
    assert(isize[0]==loc[0] && isize[1]==loc[1] && isize[2]==loc[2]);

    #ifdef _CUDA_COMP_
    rho=(Real*)malloc(isize[0]*isize[1]*isize[2]*sizeof(Real));
    cudaMalloc((void**) &rho_gpu, isize[0]*isize[1]*isize[2]*sizeof(Real));
    //cudaMalloc((void**) &phi_gpu, isize[0]*isize[1]*isize[2]*sizeof(Real));
    cudaMalloc((void**) &phi_hat, alloc_max);
    #else
    rho=(Real*)accfft_alloc(isize[0]*isize[1]*isize[2]*sizeof(Real));
    phi_hat=(myComplex*)accfft_alloc(alloc_max);
    #endif

    #ifdef _CUDA_COMP_
      #ifndef _SP_COMP_
        plan = accfft_plan_dft_3d_r2c_gpu(n, rho_gpu, (Real*)phi_hat, c_comm, ACCFFT_MEASURE);
      #else
        plan = accfft_plan_dft_3d_r2c_gpuf(n, rho_gpu, (Real*)phi_hat, c_comm, ACCFFT_MEASURE);
      #endif
    #else
    accfft_init(desired_threads);
      #ifndef _SP_COMP_
      plan = accfft_plan_dft_3d_r2c(n, rho, (Real*)phi_hat, c_comm, ACCFFT_MEASURE);
      #else
      plan = accfft_plan_dft_3d_r2cf(n, rho, (Real*)phi_hat, c_comm, ACCFFT_MEASURE);
      #endif
    #endif

    if (TStreamer::channels != 1) {
      cout << "PoissonSolverScalar_MPI(): Error: TStreamer::channels is " << TStreamer::channels << " (should be 1).\n";
      abort();
    }
  }

  void solve()
  {
    #ifdef _CUDA_COMP_
    cudaMemcpy(rho_gpu, rho, isize[0]*isize[1]*isize[2]*sizeof(Real),
              cudaMemcpyHostToDevice);
    #endif

    // Perform forward FFT
    MPI_Barrier(c_comm);

    #ifdef _CUDA_COMP_
      #ifndef _SP_COMP_
      accfft_execute_r2c_gpu(plan,rho_gpu,phi_hat);
      #else
      accfft_execute_r2c_gpuf(plan,rho_gpu,phi_hat);
      #endif
    #else
      #ifndef _SP_COMP_
      accfft_execute_r2c(plan,rho,phi_hat);
      #else
      accfft_execute_r2cf(plan,rho,phi_hat);
      #endif
    #endif

    // Spectral solve
    MPI_Barrier(c_comm);
    const double h = grid.getBlocksInfo().front().h_gridpoint;
      #ifdef _CUDA_COMP_
      const int NN[3] = {totbpd[0]*bs[0], totbpd[1]*bs[1], totbpd[2]*bs[2]};
      _fourier_filter_gpu(phi_hat, NN, osize, ostart, h);
      #else
      _fourier_filter(phi_hat, h);
    #endif

    // Spectral solve
    MPI_Barrier(c_comm);
    // Perform backward FFT
    #ifdef _CUDA_COMP_
      #ifndef _SP_COMP_
      accfft_execute_c2r_gpu(plan,phi_hat,rho_gpu);
      //accfft_execute_c2r_gpu(plan,phi_hat,phi_gpu);
      #else
      accfft_execute_c2r_gpuf(plan,phi_hat,rho_gpu);
      //accfft_execute_c2r_gpuf(plan,phi_hat,phi_gpu);
      #endif
    #else
      #ifndef _SP_COMP_
      accfft_execute_c2r(plan,phi_hat,rho);
      #else
      accfft_execute_c2rf(plan,phi_hat,rho);
      #endif
    #endif

    #ifdef _CUDA_COMP_
    cudaMemcpy(rho, rho_gpu, isize[0]*isize[1]*isize[2]*sizeof(Real),
            cudaMemcpyDeviceToHost);
    #endif
    _fft2cub(rho);
  }

  ~PoissonSolverScalarFFTW_ACC()
  {
    #ifndef _CUDA_COMP_
    accfft_free(rho);
    accfft_free(phi_hat);
    accfft_destroy_plan(plan);
    accfft_cleanup();
    #else
    free(rho);
    cudaFree(rho_gpu);
    cudaFree(phi_hat);
    accfft_destroy_plan_gpu(plan);
      #ifndef _SP_COMP_
    accfft_cleanup_gpu();
    #else
    accfft_cleanup_gpuf();
    #endif
    #endif
    MPI_Comm_free(&c_comm);
  }

  inline size_t _offset(const int blockID) const
  {
    const BlockInfo &info = local_infos[blockID];
    return _offset(info);
  }
  inline size_t _offset_ext(const BlockInfo &info) const
  {
    assert(local_infos[info.blockID].blockID == info.blockID);
    return _offset(local_infos[info.blockID]);
    //for(const auto & local_info : local_infos)
    //  if(local_info.blockID == info.blockID)
    //    return _offset(local_info);
    //printf("PSolver cannot find obstacle block\n");
    //abort();
  }
  inline size_t _offset(const BlockInfo &info) const
  {
    const int myIstart[3] = {
      info.index[0]*bs[0],
      info.index[1]*bs[1],
      info.index[2]*bs[2]
    };
    return myIstart[2] + myN[2]*myIstart[1] + myN[2]*myN[1]*myIstart[0];
  }
  inline size_t _dest(const size_t offset,const int z,const int y,const int x) const
  {
    return offset + z + myN[2] * (y + myN[1] * x);
  }
  inline void _cub2fftw(const size_t offset, const int z, const int y, const int x, const Real rhs) const
  {
    const size_t dest_index = _dest(offset, z, y, x);
    assert(dest_index >= 0 && dest_index < myN[0] * myN[1] * myN[2]);
    rho[dest_index] = rhs;
  }
};
