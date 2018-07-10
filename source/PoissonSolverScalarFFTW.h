//
//  CubismUP_3D
//
//  Written by Guido Novati ( novatig@ethz.ch ).
//  Copyright (c) 2017 ETHZ. All rights reserved.
//

#pragma once

#include "Definitions.h"

#include <fftw3.h>
#include <fftw3-mpi.h>

#ifndef _FLOAT_PRECISION_
typedef fftw_complex mycomplex;
typedef fftw_plan myplan;
#else
typedef fftwf_complex mycomplex;
typedef fftwf_plan myplan;
#endif

using namespace std;

#include <BlockInfo.h>

template<typename TGrid, typename TStreamer>
class PoissonSolverScalarFFTW_MPI
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
  const size_t nz_hat = gsize[2]/2+1;
  const double norm_factor = 1./(gsize[0]*gsize[1]*gsize[2]);
  const double h = grid.getBlocksInfo().front().h_gridpoint;
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
        assert(src_index>=0 && src_index<gsize[0]*gsize[1]*nz_hat*2);
        b(ix,iy,iz).p = out[src_index];
      }
    }
  }

  void _solve()
  {
    mycomplex *const in_out = (mycomplex *) data;
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
        const Real kinv = kx || ky || kz ? -1/(rkx*rkx+rky*rky+rkz*rkz) : 0;
        in_out[linidx][0] *= kinv*norm_factor;
        in_out[linidx][1] *= kinv*norm_factor;
      }
    #endif

    //this is sparta!
    if (local_1_start == 0) in_out[0][0] = in_out[0][1] = 0;
  }

 public:
  Real * data = nullptr;

  PoissonSolverScalarFFTW_MPI(TGrid& g): grid(g)
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
    #ifndef _FLOAT_PRECISION_
      const int retval = fftw_init_threads();
    #else
      const int retval = fftwf_init_threads();
    #endif
    if(retval==0) {
      cout << "FFTWBase::setup(): Call to fftw_init_threads() returned zero.\n";
      abort();
    }

    const int desired_threads = omp_get_max_threads();
    #ifndef _FLOAT_PRECISION_
      fftw_plan_with_nthreads(desired_threads);
      fftw_mpi_init();
      alloc_local = fftw_mpi_local_size_3d_transposed(
        gsize[0], gsize[1], gsize[2]/2+1, comm,
        &local_n0, &local_0_start, &local_n1, &local_1_start);

      data = fftw_alloc_real(2*alloc_local);
      fwd = fftw_mpi_plan_dft_r2c_3d(gsize[0], gsize[1], gsize[2],
        data, (mycomplex *)data, comm, FFTW_MPI_TRANSPOSED_OUT | FFTW_MEASURE);
      bwd = fftw_mpi_plan_dft_c2r_3d(gsize[0], gsize[1], gsize[2],
        (mycomplex *)data, data, comm, FFTW_MPI_TRANSPOSED_IN  | FFTW_MEASURE);
    #else
      fftwf_plan_with_nthreads(desired_threads);
      fftwf_mpi_init();
      alloc_local = fftwf_mpi_local_size_3d_transposed(
        gsize[0], gsize[1], gsize[2]/2+1, comm,
        &local_n0, &local_0_start, &local_n1, &local_1_start);

      data = fftwf_alloc_real(2*alloc_local);
      fwd = fftwf_mpi_plan_dft_r2c_3d(gsize[0], gsize[1], gsize[2],
        data, (mycomplex *)data, comm, FFTW_MPI_TRANSPOSED_OUT | FFTW_MEASURE);
      bwd = fftwf_mpi_plan_dft_c2r_3d(gsize[0], gsize[1], gsize[2],
        (mycomplex *)data, data, comm, FFTW_MPI_TRANSPOSED_IN  | FFTW_MEASURE);
    #endif
    //std::cout <<    bs[0] << " " <<    bs[1] << " " <<    bs[2] << " ";
    //std::cout <<   myN[0] << " " <<   myN[1] << " " <<   myN[2] << " ";
    //std::cout << gsize[0] << " " << gsize[1] << " " << gsize[2] << " ";
    //std::cout << mybpd[0] << " " << mybpd[1] << " " << mybpd[2] << std::endl;
  }

  void solve()
  {
    //_cub2fftw(data);

    #ifndef _FLOAT_PRECISION_
    fftw_execute(fwd);
    #else
    fftwf_execute(fwd);
    #endif

    _solve();

    #ifndef _FLOAT_PRECISION_
    fftw_execute(bwd);
    #else
    fftwf_execute(bwd);
    #endif

    _fftw2cub(data);
  }

  void dispose()
  {
    #ifndef _FLOAT_PRECISION_
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
    return myIstart[2] +2*nz_hat*(myIstart[1]+ myN[1]*myIstart[0]);
  }
  inline size_t _dest(const size_t offset,const int z,const int y,const int x) const
  {
    return offset + z + 2*nz_hat*(y + myN[1] * x);
  }
  inline void _cub2fftw(const size_t offset, const int z, const int y, const int x, const Real rhs) const
  {
    const size_t dest_index = _dest(offset, z, y, x);
    assert(dest_index>=0 && dest_index<gsize[0]*gsize[1]*nz_hat*2);
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
        assert(dest_index>=0 && dest_index<gsize[0]*gsize[1]*nz_hat*2);
        out[dest_index] = b(ix,iy,iz).p;
        //TStreamer::operate(b.data[iz][iy][ix], &out[dest_index]);
      }
    }
  }
};
