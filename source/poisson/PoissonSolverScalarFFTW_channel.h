       //
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#pragma once

#include "../Definitions.h"

#include <fftw3.h>
#include <fftw3-mpi.h>

#ifndef CUP_SINGLE_PRECISION
typedef fftw_complex mycomplex;
typedef fftw_plan myplan;
#else
typedef fftwf_complex mycomplex;
typedef fftwf_plan myplan;
#endif

using namespace std;

#include <Cubism/BlockInfo.h>

template<typename TGrid, typename TStreamer>
class PoissonSolverScalarFFTW_channel
{
  typedef typename TGrid::BlockType BlockType;
  TGrid& grid;
  //mycomplex local_rhs, local_work;
  myplan fwd, bwd;

  const int bs[3] = {BlockType::sizeX, BlockType::sizeY, BlockType::sizeZ};
  const vector<BlockInfo> local_infos = grid.getResidentBlocksInfo();
  const size_t mybpd[3] = {
      static_cast<size_t>(grid.getResidentBlocksPerDimension(0)),
      static_cast<size_t>(grid.getResidentBlocksPerDimension(1)),
      static_cast<size_t>(grid.getResidentBlocksPerDimension(2))
  };
  const size_t gsize[3] = {
      static_cast<size_t>(grid.getBlocksPerDimension(0)*bs[0]),
      static_cast<size_t>(grid.getBlocksPerDimension(1)*bs[1]),
      static_cast<size_t>(grid.getBlocksPerDimension(2)*bs[2])
  };
  const size_t myN[3]={ mybpd[0]*bs[0], mybpd[1]*bs[1], mybpd[2]*bs[2] };
  const double h = grid.getBlocksInfo().front().h_gridpoint;
  const double norm_factor = .125/(gsize[0]*h*gsize[1]*h*gsize[2]*h);
  ptrdiff_t alloc_local=0, local_n0=0, local_0_start=0, local_n1=0, local_1_start=0;

 protected:

  void _fftw2cub(const Real * const out) const
  {
    #pragma omp parallel for
    for(size_t i=0; i<local_infos.size(); ++i) {
      const BlockInfo info = local_infos[i];
      BlockType& b = *(BlockType*)info.ptrBlock;
      const size_t offset = _offset(info);

      for(int ix=0; ix<BlockType::sizeX; ix++)
      for(int iy=0; iy<BlockType::sizeY; iy++)
      for(int iz=0; iz<BlockType::sizeZ; iz++) {
        const size_t src_index = _dest(offset, iz, iy, ix);
        assert(src_index>=0 && src_index<gsize[0]*gsize[1]*gsize[2]);
        b(ix,iy,iz).p = out[src_index];
      }
    }
  }

  void _solve()
  {
    Real *const in_out = data;
    const Real waveFactX = M_PI/(gsize[0]*h);
    const Real waveFactY = M_PI/(gsize[1]*h);
    const Real waveFactZ = M_PI/(gsize[2]*h);
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
      const long kx = i;
      const long ky = shifty + j; //memory index plus shift due to decomp
      const long kz = k;
      const Real rkx = (kx+0.5)*waveFactX;
      const Real rky = (ky+0.5)*waveFactY;
      const Real rkz = (kz+0.5)*waveFactZ;
      const Real kinv = -1/(rkx*rkx + rky*rky + rkz*rkz);
      in_out[linidx] *= kinv*norm_factor;
    }
    if (shifty == 0) in_out[0] = 0;
  }

 public:
  Real * data = nullptr;

  PoissonSolverScalarFFTW_channel(TGrid& g): grid(g)
  {
    if (TStreamer::channels != 1) {
      cout << "PoissonSolverScalar_MPI(): Error: TStreamer::channels is "
           << TStreamer::channels << " (should be 1).\n";
      abort();
    }
    MPI_Comm comm = grid.getCartComm();

    {
      int supported_threads;
      MPI_Query_thread(&supported_threads);
      if (supported_threads<MPI_THREAD_FUNNELED) {
        cout << "MPI implementation does not support threads.\n";
        abort();
      }
    }
    #ifndef CUP_SINGLE_PRECISION
      const int retval = fftw_init_threads();
    #else
      const int retval = fftwf_init_threads();
    #endif
    if(retval==0) {
      cout << "FFTWBase::setup(): Call to fftw_init_threads() returned zero.\n";
      abort();
    }

    const int desired_threads = omp_get_max_threads();
    #ifndef CUP_SINGLE_PRECISION
      fftw_plan_with_nthreads(desired_threads);
      fftw_mpi_init();

      alloc_local = fftw_mpi_local_size_3d_transposed(
        gsize[0], gsize[1], gsize[2], comm,
        &local_n0, &local_0_start, &local_n1, &local_1_start);
      data = fftw_alloc_real(alloc_local);
      fwd = fftw_mpi_plan_r2r_3d(gsize[0], gsize[1], gsize[2],
        data, data, comm, FFTW_REDFT10, FFTW_REDFT10, FFTW_REDFT10, FFTW_MPI_TRANSPOSED_OUT  | FFTW_MEASURE);
      bwd = fftw_mpi_plan_r2r_3d(gsize[0], gsize[1], gsize[2],
        data, data, comm, FFTW_REDFT01, FFTW_REDFT01, FFTW_REDFT01, FFTW_MPI_TRANSPOSED_IN  | FFTW_MEASURE);
    #else
      fftwf_plan_with_nthreads(desired_threads);
      fftwf_mpi_init();
      alloc_local = fftwf_mpi_local_size_3d_transposed(
        gsize[0], gsize[1], gsize[2], comm,
        &local_n0, &local_0_start, &local_n1, &local_1_start);
      data = fftwf_alloc_real(alloc_local);
      fwd = fftwf_mpi_plan_r2r_3d(gsize[0], gsize[1], gsize[2],
        data, data, comm, FFTW_REDFT10, FFTW_REDFT10, FFTW_REDFT10, FFTW_MPI_TRANSPOSED_OUT  | FFTW_MEASURE);
      bwd = fftwf_mpi_plan_r2r_3d(gsize[0], gsize[1], gsize[2],
        data, data, comm, FFTW_REDFT01, FFTW_REDFT01, FFTW_REDFT01, FFTW_MPI_TRANSPOSED_IN  | FFTW_MEASURE);
    #endif
    //std::cout <<    bs[0] << " " <<    bs[1] << " " <<    bs[2] << " ";
    //std::cout <<   myN[0] << " " <<   myN[1] << " " <<   myN[2] << " ";
    //std::cout << gsize[0] << " " << gsize[1] << " " << gsize[2] << " ";
    //std::cout << mybpd[0] << " " << mybpd[1] << " " << mybpd[2] << std::endl;
  }

  void solve()
  {
    //_cub2fftw(data);
    #ifndef CUP_SINGLE_PRECISION
    fftw_execute(fwd);
    #else
    fftwf_execute(fwd);
    #endif

    _solve();

    #ifndef CUP_SINGLE_PRECISION
    fftw_execute(bwd);
    #else
    fftwf_execute(bwd);
    #endif
    _fftw2cub(data);
  }

  void dispose()
  {
    #ifndef CUP_SINGLE_PRECISION
      fftw_destroy_plan(fwd);
      fftw_destroy_plan(bwd);
      fftw_free(data);
      fftw_mpi_cleanup();
    #else
      fftwf_destroy_plan(fwd);
      fftwf_destroy_plan(bwd);
      fftwf_free(data);
      fftwf_mpi_cleanup();
    #endif
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
    return myIstart[2] +myN[2]*(myIstart[1]+ myN[1]*myIstart[0]);
  }
  inline size_t _dest(const size_t offset,const int z,const int y,const int x) const
  {
    return offset + z + myN[2]*(y + myN[1] * x);
  }
  inline void _cub2fftw(const size_t offset, const int z, const int y, const int x, const Real rhs) const
  {
    const size_t dest_index = _dest(offset, z, y, x);
    assert(dest_index>=0 && dest_index<gsize[0]*gsize[1]*myN[2]);
    data[dest_index] = rhs;
  }

  void _cub2fftw(Real * const out) const
  {
    #pragma omp parallel for
    for(size_t i=0; i<local_infos.size(); ++i) {
      const BlockInfo info = local_infos[i];
      BlockType& b = *(BlockType*)info.ptrBlock;
      const size_t offset = _offset(info);

      for(int ix=0; ix<BlockType::sizeX; ix++)
      for(int iy=0; iy<BlockType::sizeY; iy++)
      for(int iz=0; iz<BlockType::sizeZ; iz++) {
        const size_t dest_index = _dest(offset, iz, iy, ix);
        assert(dest_index>=0 && dest_index<gsize[0]*gsize[1]*myN[2]);
        out[dest_index] = b(ix,iy,iz).p;
        //TStreamer::operate(b.data[iz][iy][ix], &out[dest_index]);
      }
    }
  }
};
