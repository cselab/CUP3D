//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch) and Christian Conti.
//

#ifndef CubismUP_3D_ProcessOperators_h
#define CubismUP_3D_ProcessOperators_h

#include "SimulationData.h"

CubismUP_3D_NAMESPACE_BEGIN

// -gradp, divergence, advection
template<typename Lab, typename Kernel>
void processOMP(double dt, const SimulationData & sim)
{
  const Kernel kernel(dt);
  SynchronizerMPI<Real>& Synch = sim.grid->sync(kernel);
  const int nthreads = omp_get_max_threads();
  LabMPI * labs = new LabMPI[nthreads];
  #pragma omp parallel for schedule(static, 1)
  for(int i = 0; i < nthreads; ++i) {
    labs[i].setBC(sim.BCx_flag, sim.BCy_flag, sim.BCz_flag);
    labs[i].prepare(* sim.grid, Synch);
  }

  MPI_Barrier(sim.grid->getCartComm());
  std::vector<BlockInfo> avail0 = Synch.avail_inner();
  const int Ninner = avail0.size();

  #pragma omp parallel
  {
    int tid = omp_get_thread_num();
    LabMPI& lab = labs[tid];

    #pragma omp for schedule(static)
    for(int i=0; i<Ninner; i++) {
      BlockInfo info = avail0[i];
      FluidBlock& b = *(FluidBlock*)info.ptrBlock;
      lab.load(info, 0);
      kernel(lab, info, b); // why is this using the local blockInfo? or is it global? is dh correct?
    }
  }

  std::vector<BlockInfo> avail1 = Synch.avail_halo();
  const int Nhalo = avail1.size();

  #pragma omp parallel
  {
    int tid = omp_get_thread_num();
    LabMPI& lab = labs[tid];

    #pragma omp for schedule(static)
    for(int i=0; i<Nhalo; i++) {
        BlockInfo info = avail1[i];
        FluidBlock& b = *(FluidBlock*)info.ptrBlock;
        lab.load(info, 0);
        kernel(lab, info, b);
    }
  }

  if(labs!=NULL) {
    delete [] labs;
    labs=NULL;
  }

  MPI_Barrier(sim.grid->getCartComm());
}

inline Real findMaxUOMP(
        const std::vector<BlockInfo> myInfo,
        FluidGridMPI& grid,
        const Real*const uinf
      )
{
  Real maxU = 0;
  const int N = myInfo.size();
  const Real U[3] = {uinf[0], uinf[1], uinf[2]};
  #pragma omp parallel for schedule(static) reduction(max : maxU)
  for(int i=0; i<N; i++)
  {
    const BlockInfo& info = myInfo[i];
    const FluidBlock& b = *(const FluidBlock *)info.ptrBlock;

    for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
    for(int iy=0; iy<FluidBlock::sizeY; ++iy)
    for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
      Real u=b(ix,iy,iz).u+U[0], v=b(ix,iy,iz).v+U[1], w=b(ix,iy,iz).w+U[2];
      Real au = std::fabs(u), av = std::fabs(v), aw = std::fabs(w);
      const Real maxUl = std::max({ au, av, aw});
      maxU = std::max(maxU, maxUl);
    }
  }
  return maxU;
}

#ifdef CUP_ASYNC_DUMP
static void copyDumpGrid(FluidGridMPI& grid, DumpGridMPI& dump)
{
  std::vector<BlockInfo> vInfo1 = grid.getBlocksInfo();
  std::vector<BlockInfo> vInfo2 = dump.getBlocksInfo();
  const int N = vInfo1.size();
  if(vInfo1.size() != vInfo2.size()) {
     printf("Async dump fail 1.\n");
     fflush(0);
     MPI_Abort(grid.getCartComm(), MPI_ERR_OTHER);
   }
  #pragma omp parallel for schedule(static)
  for(int i=0; i<N; i++) {
    const BlockInfo& info1 = vInfo1[i];
    const BlockInfo& info2 = vInfo2[i];

    #ifndef NDEBUG
      Real p1[3], p2[3];
      info1.pos(p1, 0,0,0);
      info2.pos(p2, 0,0,0);
      if (fabs(p1[0]-p2[0])>info1.h_gridpoint/2 ||
          fabs(p1[1]-p2[1])>info1.h_gridpoint/2 ||
          fabs(p1[2]-p2[2])>info1.h_gridpoint/2) {
             printf("Async dump fail 2.\n");
             fflush(0);
             MPI_Abort(grid.getCartComm(), MPI_ERR_OTHER);
          }
    #endif

    const FluidBlock& b = *(FluidBlock*)info1.ptrBlock;
     DumpBlock& d = *( DumpBlock*)info2.ptrBlock;
    for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
    for(int iy=0; iy<FluidBlock::sizeY; ++iy)
    for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
      d(ix,iy,iz).u = b(ix,iy,iz).u;
      d(ix,iy,iz).v = b(ix,iy,iz).v;
      d(ix,iy,iz).w = b(ix,iy,iz).w;
      d(ix,iy,iz).chi = b(ix,iy,iz).chi;
      d(ix,iy,iz).p = b(ix,iy,iz).p;
    }
  }
}
#endif

CubismUP_3D_NAMESPACE_END
#endif // CubismUP_3D_ProcessOperators_h
