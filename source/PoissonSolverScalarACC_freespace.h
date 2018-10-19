//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#pragma once

#include "Definitions.h"
typedef Real myCmpl[2];

/******************************************************************************/
void initGreen(const int nx, const int ny, const int nz, const Real h, Real*const m_kernel, const MPI_Comm comm);
void dSolveFreespace(void*const P, const int nx,const int ny,const int nz,
  const int locx,const int locy,const int locz, const MPI_Comm comm,
  const int ox,const int oy,const int oz, const Real*const G_hat,
  Real*const cub_rhs, Real*const fft_rhs, Real*const gpu_rhs, myCmpl*const phi);
MPI_Comm my_accfft_create_comm(MPI_Comm C, int c_dims[2]);
Real* my_cudaMalloc(const size_t size);
void my_cudaFree(Real* const ptr);
void my_accfft_cleanup(void* const plan);
size_t my_accfft_local_size(int M[3], int isz[3], int ist[3], int osz[3],
  int ost[3], MPI_Comm c_comm);
void* my_accfft_plan_dft(int M[3], Real* gpurhs, MPI_Comm c_comm);
void my_cudaMemset_zero(Real* const ptr, const size_t size);

template<typename TGrid, typename TStreamer>
class PoissonSolverScalarFFTW_ACC
{
  typedef typename TGrid::BlockType BlockType;
  TGrid& grid;

  const int bs[3] = {BlockType::sizeX, BlockType::sizeY, BlockType::sizeZ};
  const int mybpd[3] = {  grid.getResidentBlocksPerDimension(0),
                          grid.getResidentBlocksPerDimension(1),
                          grid.getResidentBlocksPerDimension(2) };
  const int totbpd[3] = { grid.getBlocksPerDimension(0),
                          grid.getBlocksPerDimension(1),
                          grid.getBlocksPerDimension(2) };

  const size_t myN[3] = {
    static_cast<size_t>(mybpd[0]*bs[0]),
    static_cast<size_t>(mybpd[1]*bs[1]),
    static_cast<size_t>(mybpd[2]*bs[2])
  };
  const int totN[3] = { totbpd[0]*bs[0], totbpd[1]*bs[1], totbpd[2]*bs[2] };
  const int nprocs = getSize( grid.getCartComm() );
  const int procid = getRank( grid.getCartComm() );

  int isz[3]={0,0,0}, osz[3]={0,0,0}, ist[3]={0,0,0}, ost[3]={0,0,0};

  const vector<BlockInfo> local_infos = grid.getResidentBlocksInfo();
  const size_t N = local_infos.size();

  const int mx = 2 * totN[0] - 1, my = 2 * totN[1] - 1, mz = 2 * totN[2] - 1;
  const int mz_pad = (mz/2 +1)*2, myftNx = (mx+1)/nprocs;

  int c_dims[2] = { nprocs, 1 };
  MPI_Comm c_comm = my_accfft_create_comm( grid.getCartComm(), c_dims);

  int M[3] = {mx, my, mz};
  const size_t alloc_max = my_accfft_local_size(M, isz,ist,osz,ost, c_comm);
  Real * const cub_rhs = (Real*) malloc(myN[0]* myN[1]* myN[2] * sizeof(Real) );
  Real * const fft_rhs = (Real*) malloc(myftNx*totN[1]*totN[2] * sizeof(Real) );
  Real * const gpuGhat = my_cudaMalloc(alloc_max / 2);
  Real * const gpu_rhs = my_cudaMalloc( isz[0]*isz[1]*isz[2]*sizeof(Real) );
  myCmpl * const gpu_phi = my_cudaMalloc(alloc_max);
  void * const plan = my_accfft_plan_dft(M, gpu_rhs, c_comm);

  static int getSize(MPI_Comm comm) {
    int ret; MPI_Comm_size(comm, &ret); return ret;
  }
  static int getRank(MPI_Comm comm) {
    int ret; MPI_Comm_rank(comm, &ret); return ret;
  }
public:
  PoissonSolverScalarFFTW_ACC(TGrid& g) : grid(g) {
   printf("[mpi rank %d] isize  %3d %3d %3d osize  %3d %3d %3d\n", procid,
        isz[0],isz[1],isz[2], osz[0],osz[1],osz[2] );
   printf("[mpi rank %d] istart %3d %3d %3d ostart %3d %3d %3d\n", procid,
        ist[0],ist[1],ist[2], ost[0],ost[1],ost[2] );
   initGreen(totN[0],totN[1],totN[2],local_infos[0].h_gridpoint,gpuGhat,c_comm);
   if(procid<nprocs/2) {
     assert(isz[0]==myftNx);
     assert(isz[1]==totN[1]);
     assert(isz[2]==totN[2]);
   }
  }

  void solve() {
    const MPI_Comm cart_comm = grid.getCartComm();
    memset(fft_rhs, 0, myftNx*totN[1]*totN[2] * sizeof(Real) );
    my_cudaMemset_zero(gpu_rhs, alloc_max);
    dSolveFreespace(plan, totN[0],totN[1],totN[2], myN[0],myN[1],myN[2],
     cart_comm, osz[0],osz[1],osz[2], gpuGhat,cub_rhs,fft_rhs, gpu_rhs,gpu_phi);

    _fft2cub();
  }

  ~PoissonSolverScalarFFTW_ACC() {
    free(cub_rhs);
    free(fft_rhs);
    my_cudaFree(gpu_rhs);
    my_cudaFree(gpuGhat);
    my_accfft_cleanup(plan);
    MPI_Comm_free(&c_comm);
  }

  inline size_t _offset(const int blockID) const {
    const BlockInfo &info = local_infos[blockID];
    return _offset(info);
  }
  inline size_t _offset_ext(const BlockInfo &info) const {
    assert(local_infos[info.blockID].blockID == info.blockID);
    return _offset(local_infos[info.blockID]);
    //for(const auto & local_info : local_infos)
    //  if(local_info.blockID == info.blockID)
    //    return _offset(local_info);
    //printf("PSolver cannot find obstacle block\n");
    //abort();
  }
  inline size_t _offset(const BlockInfo &info) const {
    const int myIstart[3] = {
      info.index[0]*bs[0],
      info.index[1]*bs[1],
      info.index[2]*bs[2]
    };
    return myIstart[2] + myN[2]*myIstart[1] + myN[2]*myN[1]*myIstart[0];
  }
  inline size_t _dest(const size_t offset,const int z,const int y,const int x) const {
    return offset + z + myN[2] * (y + myN[1] * x);
  }
  inline void _cub2fftw(const size_t offset, const int z, const int y, const int x, const Real rhs) const {
    const size_t dest_index = _dest(offset, z, y, x);
    assert(dest_index >= 0 && dest_index < myN[0] * myN[1] * myN[2]);
    cub_rhs[dest_index] = rhs;
  }
  void _fft2cub() const {
    #pragma omp parallel for schedule(static)
    for(size_t i=0; i<N; ++i) {
      const BlockInfo info = local_infos[i];
      BlockType& b = *(BlockType*)info.ptrBlock;
      const size_t offset = _offset(info);
      for(int ix=0; ix<bs[0]; ix++)
      for(int iy=0; iy<bs[1]; iy++)
      for(int iz=0; iz<bs[2]; iz++) {
        const size_t src_index = _dest(offset, iz, iy, ix);
        assert(src_index >= 0 && src_index < myN[0] * myN[1] * myN[2]);
        b(ix,iy,iz).p = cub_rhs[src_index];
      }
    }
  }
};
