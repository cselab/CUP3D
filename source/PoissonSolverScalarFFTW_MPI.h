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

#include <fftw3.h>
#include <fftw3-mpi.h>

#ifndef _SP_COMP_
typedef fftw_complex mycomplex;
typedef fftw_plan myplan;
#else
typedef fftwf_complex mycomplex;
typedef fftwf_plan myplan;
#endif

using namespace std;

#include <BlockInfo.h>

class FFTWBase_MPI
{
  static int registered_objects;
  static bool initialized; //a la singleton

  static void _setup(const int desired_threads)
  {
    if (!initialized) {
      initialized = true;
      int supported_threads;
      MPI_Query_thread(&supported_threads);
      if (supported_threads>=MPI_THREAD_FUNNELED) {
        #ifndef _SP_COMP_
          const int retval = fftw_init_threads();
        #else
          const int retval = fftwf_init_threads();
        #endif
        if(retval==0) {
          cout << "FFTWBase::setup(): Oops the call to fftw_init_threads() returned zero. Aborting\n";
          abort();
        }
        else {
          #ifndef _SP_COMP_
            fftw_plan_with_nthreads(desired_threads);
          #else
            fftwf_plan_with_nthreads(desired_threads);
          #endif
        }
      }

      #ifndef _SP_COMP_
      fftw_mpi_init();
      #else
      fftwf_mpi_init();
      #endif
    }

    registered_objects++;
  }

 public:

  FFTWBase_MPI(const int desired_threads) { _setup(desired_threads); }

  static void dispose()
  {
    registered_objects--;

    if (registered_objects == 0) {
      #ifndef _SP_COMP_
      fftw_mpi_cleanup();
      #else
      fftwf_mpi_cleanup();
      #endif
    }
  }
};

template<typename TGrid, typename TStreamer>
class PoissonSolverScalarFFTW_MPI : FFTWBase_MPI
{
  typedef typename TGrid::BlockType BlockType;
  TGrid& grid;

  bool initialized;
  //mycomplex local_rhs, local_work;
  myplan fwd, bwd;

  const int bs[3] = {BlockType::sizeX, BlockType::sizeY, BlockType::sizeZ};
  const vector<BlockInfo> local_infos = grid.getResidentBlocksInfo();
  const size_t N = local_infos.size();
  const size_t mybpd[3] = {
      static_cast<size_t>(grid.getResidentBlocksPerDimension(0)),
      static_cast<size_t>(grid.getResidentBlocksPerDimension(1)),
      static_cast<size_t>(grid.getResidentBlocksPerDimension(2))
  };
  const size_t gsize[3] = {
      grid.getBlocksPerDimension(0)*bs[0],
      grid.getBlocksPerDimension(1)*bs[1],
      grid.getBlocksPerDimension(2)*bs[2]
  };
  const size_t myN[3]={ mybpd[0]*bs[0], mybpd[1]*bs[1], mybpd[2]*bs[2] };
  const size_t nz_hat = gsize[2]/2+1;
  ptrdiff_t alloc_local=0, local_n0=0, local_0_start=0, local_n1=0, local_1_start=0;

protected:

  void _setup(Real*& rhs,const size_t nx,const size_t ny,const size_t nz,MPI_Comm comm)
  {
    if (!initialized) {
      initialized = true;
      #ifndef _SP_COMP_
      alloc_local = fftw_mpi_local_size_3d_transposed(nx, ny, nz/2+1, comm,
                        &local_n0, &local_0_start, &local_n1, &local_1_start);
      rhs = fftw_alloc_real(2*alloc_local);
      fwd = fftw_mpi_plan_dft_r2c_3d(nx, ny, nz, rhs, (mycomplex *)rhs, comm,
                                      FFTW_MPI_TRANSPOSED_OUT | FFTW_MEASURE);
      bwd = fftw_mpi_plan_dft_c2r_3d(nx, ny, nz, (mycomplex *)rhs, rhs, comm,
                                      FFTW_MPI_TRANSPOSED_IN  | FFTW_MEASURE);
      #else
      alloc_local = fftwf_mpi_local_size_3d_transposed(nx, ny, nz/2+1, comm,
                          &local_n0, &local_0_start, &local_n1, &local_1_start);
      rhs = fftwf_alloc_real(2*alloc_local);
      fwd = fftwf_mpi_plan_dft_r2c_3d(nx, ny, nz, rhs, (mycomplex *)rhs, comm,
                                      FFTW_MPI_TRANSPOSED_OUT | FFTW_MEASURE);
      bwd = fftwf_mpi_plan_dft_c2r_3d(nx, ny, nz, (mycomplex *)rhs, rhs, comm,
                                      FFTW_MPI_TRANSPOSED_IN  | FFTW_MEASURE);
      #endif
    }
  }

  void _fftw2cub(const Real * const out) const
  {
    #pragma omp parallel for
    for(int i=0; i<N; ++i) {
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

  void _solve(mycomplex* in_out,const size_t nx,const size_t ny, const size_t
    nz, const size_t _nz_hat, const double norm_factor, const double h)
  {
    #if 0
    const Real h2 = h*h;
    const Real factor = h2*norm_factor;

    #pragma omp parallel for
    for(int j=0; j<local_n1; ++j)
    for(int i=0; i<nx; ++i)
    for(int k = 0; k<nz_hat; ++k) {
      const int linidx = (j*nx+i)*nz_hat + k;
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
    const Real waveFactX = 2.0*M_PI/(nx*h);
    const Real waveFactY = 2.0*M_PI/(ny*h);
    const Real waveFactZ = 2.0*M_PI/(nz*h);

    #pragma omp parallel for
    for(int j = 0; j<local_n1; ++j)
    for(int i=0; i<nx; ++i)
    for(int k = 0; k<nz_hat; ++k) {
      const int linidx = (j*nx+i)*nz_hat + k;
      //assert(linidx >=0 && linidx<nx*ny*nz_hat); // linking error with openmp
      const int kx=(i <= nx/2) ? i : -(nx-i);
      const int ky=(local_1_start+j<=ny/2)?local_1_start+j : local_1_start+j-ny;
      const int kz=(k <= nz/2) ? k : -(nz-k);
      const Real rkx = kx*waveFactX, rky = ky*waveFactY, rkz = kz*waveFactZ;
      const Real kinv = kx || ky || kz ? -1/Real(rkx*rkx+rky*rky+rkz*rkz) : 0;
      in_out[linidx][0] *= kinv*norm_factor;
      in_out[linidx][1] *= kinv*norm_factor;
    }
    #endif

    //this is sparta!
    if (local_1_start == 0) in_out[0][0] = in_out[0][1] = 0;
  }

public:
  Real * data = nullptr;

  PoissonSolverScalarFFTW_MPI(const int desired_threads, TGrid& g):
  FFTWBase_MPI(desired_threads), grid(g), initialized(false)
  {
    if (TStreamer::channels != 1) {
      cout << "PoissonSolverScalar_MPI(): Error: TStreamer::channels is "
           << TStreamer::channels << " (should be 1).\n";
      abort();
    }
    MPI_Comm comm = grid.getCartComm();
    _setup(data, gsize[0], gsize[1], gsize[2],comm);
  }

  void solve()
  {
    #ifndef _SP_COMP_
    fftw_execute(fwd);
    #else
    fftwf_execute(fwd);
    #endif

    const double norm_factor = 1./(gsize[0]*gsize[1]*gsize[2]);
    const double h = grid.getBlocksInfo().front().h_gridpoint;
    _solve((mycomplex *)data, gsize[0], gsize[1], gsize[2], gsize[2]/2+1, norm_factor, h);

    #ifndef _SP_COMP_
    fftw_execute(bwd);
    #else
    fftwf_execute(bwd);
    #endif

    _fftw2cub(data);
  }

  void dispose()
  {
    if (initialized) {
      initialized = false;

      #ifndef _SP_COMP_
      fftw_destroy_plan(fwd);
      fftw_destroy_plan(bwd);
      fftw_free(data);
      #else
      fftwf_destroy_plan(fwd);
      fftwf_destroy_plan(bwd);
      fftwf_free(data);
      #endif
      FFTWBase_MPI::dispose();
    }
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
    for(int i=0; i<N; ++i) {
      const BlockInfo info = local_infos[i];
      BlockType& b = *(BlockType*)info.ptrBlock;
      const size_t offset = _offset(info);

      for(int ix=0; ix<BlockType::sizeX; ix++)
      for(int iy=0; iy<BlockType::sizeY; iy++)
      for(int iz=0; iz<BlockType::sizeZ; iz++) {
        const size_t dest_index = _dest(offset, iz, iy, ix);
        assert(dest_index>=0 && dest_index<gsize[0]*gsize[1]*nz_hat*2);
        TStreamer::operate(b.data[iz][iy][ix], &out[dest_index]);
      }
    }
  }
};
