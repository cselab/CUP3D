//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch) and Christian Conti.
//

#ifndef CubismUP_3D_ProcessOperators_h
#define CubismUP_3D_ProcessOperators_h

#include "Definitions.h"

// -gradp, divergence, advection
template<typename Lab, typename Kernel>
void processOMP(double dt, vector<BlockInfo>& vInfo, FluidGridMPI & grid)
{
  const Kernel kernel(dt);
  SynchronizerMPI& Synch = grid.sync(kernel);

  const int nthreads = omp_get_max_threads();
  LabMPI * labs = new LabMPI[nthreads];
  for(int i = 0; i < nthreads; ++i)
    labs[i].prepare(grid, Synch);

  MPI_Barrier(grid.getCartComm());
  vector<BlockInfo> avail0 = Synch.avail_inner();
  const int Ninner = avail0.size();

  #pragma omp parallel
  {
    int tid = omp_get_thread_num();
    LabMPI& lab = labs[tid];

    #pragma omp for schedule(dynamic,1)
    for(int i=0; i<Ninner; i++) {
      BlockInfo info = avail0[i];
      FluidBlock& b = *(FluidBlock*)info.ptrBlock;
      lab.load(info, 0);
      kernel(lab, info, b); // why is this using the local blockInfo? or is it global? is dh correct?
    }
  }

  vector<BlockInfo> avail1 = Synch.avail_halo();
  const int Nhalo = avail1.size();

  #pragma omp parallel
  {
    int tid = omp_get_thread_num();
    LabMPI& lab = labs[tid];

    #pragma omp for schedule(dynamic,1)
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

  MPI_Barrier(grid.getCartComm());
}

template<typename Lab, typename Kernel>
void processOMPold(double dt, vector<BlockInfo>& vInfo, FluidGridMPI & grid)
{
  const int N = vInfo.size();

  #pragma omp parallel
  {
    Lab lab;
    Kernel kernel(dt);
    lab.prepare(grid, kernel.stencil_start, kernel.stencil_end, true);

    #pragma omp for schedule(static)
    for (int i=0; i<N; i++) {
      BlockInfo info = vInfo[i];
      FluidBlock& b = *(FluidBlock*)info.ptrBlock;
      lab.load(info, 0);
      kernel(lab, info, b);
    }
  }
}

static Real findMaxUOMP(const vector<BlockInfo>& myInfo, FluidGridMPI& grid, const Real*const uInf)
{
  Real maxU = 0;
  const int N = myInfo.size();

  #pragma omp parallel for schedule(static) reduction(max:maxU)
  for(int i=0; i<N; i++)
  {
    const BlockInfo& info = myInfo[i];
    const FluidBlock& b = *(const FluidBlock *)info.ptrBlock;

    for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
    for(int iy=0; iy<FluidBlock::sizeY; ++iy)
    for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
      maxU = max(maxU,(Real)abs(b(ix,iy,iz).u +uInf[0]));
      maxU = max(maxU,(Real)abs(b(ix,iy,iz).v +uInf[1]));
      maxU = max(maxU,(Real)abs(b(ix,iy,iz).w +uInf[2]));
    }
  }

  return maxU;
}

#ifdef DUMPGRID
static void copyDumpGrid(FluidGridMPI& grid, DumpGridMPI& dump)
{
  vector<BlockInfo> vInfo1 = grid.getBlocksInfo();
  vector<BlockInfo> vInfo2 = dump.getBlocksInfo();
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

#endif
