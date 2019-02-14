//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#pragma once

#include "PoissonSolver.h"

template<typename TGrid, typename TStreamer>
class PoissonSolverMixed : public PoissonSolver
{
  typedef typename TGrid::BlockType BlockType;
  const SimulationData & sim;
  TGrid& grid = * sim.grid;
  //mycomplex local_rhs, local_work;
  myplan fwd, bwd;

  const int bs[3] = {BlockType::sizeX, BlockType::sizeY, BlockType::sizeZ};
  const std::vector<BlockInfo> local_infos = grid.getResidentBlocksInfo();
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
  ptrdiff_t alloc_local=0, local_n0=0, local_0_start=0, local_n1=0, local_1_start=0;

 protected:

  template<int BCX, int BCY, int BCZ> void _solve()
  {
    // if BC flag == 1 fourier, else cosine transform
    const Real normX = (BCX==1 ? 1.0 : 0.5) / ( gsize[0]*h );
    const Real normY = (BCY==1 ? 1.0 : 0.5) / ( gsize[1]*h );
    const Real normZ = (BCZ==1 ? 1.0 : 0.5) / ( gsize[2]*h );
    const Real waveFactX = (BCX==1 ? 2 : 1) * M_PI / ( gsize[0]*h );
    const Real waveFactY = (BCY==1 ? 2 : 1) * M_PI / ( gsize[1]*h );
    const Real waveFactZ = (BCZ==1 ? 2 : 1) * M_PI / ( gsize[2]*h );
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
      const long kx = BCX==1 ? ((i <= nKx/2) ? i : nKx-i) : i;
      const long ky = BCY==1 ? ((J <= nKy/2) ? J : nKy-J) : J;
      const long kz = BCZ==1 ? ((k <= nKz/2) ? k : nKz-k) : k;
      const Real rkx = ( kx + (BCX==1 ? 1.0 : 0.5) ) * waveFactX;
      const Real rky = ( ky + (BCY==1 ? 1.0 : 0.5) ) * waveFactY;
      const Real rkz = ( kz + (BCZ==1 ? 1.0 : 0.5) ) * waveFactZ;
      const Real kinv = -1/(rkx*rkx + rky*rky + rkz*rkz);
      in_out[linidx] *= kinv*norm_factor;
    }
    if (shifty == 0) in_out[0] = 0;
  }

 public:
  Real * data = nullptr;

  PoissonSolverMixed(SimulationData & s): sim(s)
  {
    if (TStreamer::channels != 1) {
      std::cout << "PoissonSolverScalar_MPI(): Error: TStreamer::channels is "
           << TStreamer::channels << " (should be 1).\n";
      abort();
    }
    MPI_Comm comm = grid.getCartComm();

    {
      int supported_threads;
      MPI_Query_thread(&supported_threads);
      if (supported_threads<MPI_THREAD_FUNNELED) {
        std::cout << "MPI implementation does not support threads.\n";
        abort();
      }
    }
    const int retval = _FFTW_(init_threads)();
    if(retval==0) {
      std::cout << "FFTWBase::setup(): Call to fftw_init_threads() returned zero." << std::endl;
      abort();
    }
    const int desired_threads = omp_get_max_threads();

    _FFTW_(plan_with_nthreads)(desired_threads);
    _FFTW_(mpi_init)();

    alloc_local = _FFTW_(mpi_local_size_3d_transposed) (
      gsize[0], gsize[1], gsize[2], comm,
      &local_n0, &local_0_start, &local_n1, &local_1_start);

    auto XplanF = sim.BCx_flag==1? FFTW_R2HC : FFTW_REDFT10;
    auto XplanB = sim.BCx_flag==1? FFTW_HC2R : FFTW_REDFT01;
    auto YplanF = sim.BCy_flag==1? FFTW_R2HC : FFTW_REDFT10;
    auto YplanB = sim.BCy_flag==1? FFTW_HC2R : FFTW_REDFT01;
    auto ZplanF = sim.BCz_flag==1? FFTW_R2HC : FFTW_REDFT10;
    auto ZplanB = sim.BCz_flag==1? FFTW_HC2R : FFTW_REDFT01;
    data = _FFTW_(alloc_real)(alloc_local);
    fwd = _FFTW_(mpi_plan_r2r_3d)(gsize[0], gsize[1], gsize[2], data, data,
      comm, XplanF, YplanF, ZplanF, FFTW_MPI_TRANSPOSED_OUT  | FFTW_MEASURE);
    bwd = _FFTW_(mpi_plan_r2r_3d)(gsize[0], gsize[1], gsize[2], data, data,
      comm, XplanB, YplanB, ZplanB, FFTW_MPI_TRANSPOSED_IN  | FFTW_MEASURE);

    //std::cout <<    bs[0] << " " <<    bs[1] << " " <<    bs[2] << " ";
    //std::cout <<   myN[0] << " " <<   myN[1] << " " <<   myN[2] << " ";
    //std::cout << gsize[0] << " " << gsize[1] << " " << gsize[2] << " ";
    //std::cout << mybpd[0] << " " << mybpd[1] << " " << mybpd[2] << std::endl;
  }

  void solve()
  {
    _cub2fftw(data);

    _FFTW_(execute)(fwd);

    if(sim.BCx_flag==1 && sim.BCy_flag==1 && sim.BCz_flag==1) _solve<1,1,1>();
    else
    if(sim.BCx_flag==1 && sim.BCy_flag==1 && sim.BCz_flag!=1) _solve<1,1,0>();
    else
    if(sim.BCx_flag==1 && sim.BCy_flag!=1 && sim.BCz_flag==1) _solve<1,0,1>();
    else
    if(sim.BCx_flag==1 && sim.BCy_flag!=1 && sim.BCz_flag!=1) _solve<1,0,0>();
    else
    if(sim.BCx_flag!=1 && sim.BCy_flag==1 && sim.BCz_flag==1) _solve<0,1,1>();
    else
    if(sim.BCx_flag!=1 && sim.BCy_flag==1 && sim.BCz_flag!=1) _solve<0,1,0>();
    else
    if(sim.BCx_flag!=1 && sim.BCy_flag!=1 && sim.BCz_flag==1) _solve<0,0,1>();
    else
    if(sim.BCx_flag!=1 && sim.BCy_flag!=1 && sim.BCz_flag!=1) _solve<0,0,0>();
    else {
      printf("Boundary conditions not recognized\n");
      abort();
    }

    _FFTW_(execute)(bwd);

    _fftw2cub(data);
  }

  ~PoissonSolverMixed()
  {
    _FFTW_(destroy_plan)(fwd);
    _FFTW_(destroy_plan)(bwd);
    _FFTW_(free)(data);
    _FFTW_(mpi_cleanup)();
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
};
