//
//  CubismUP_3D
//
//  Written by Guido Novati ( novatig@ethz.ch ).
//  This file started as an extension of code written by Christian Conti
//  Copyright (c) 2017 ETHZ. All rights reserved.
//

#ifndef CubismUP_3D_GenericCoordinator_h
#define CubismUP_3D_GenericCoordinator_h

#include "Definitions.h"

class GenericCoordinator
{
protected:
  FluidGridMPI * grid;
  vector<BlockInfo> vInfo;

  inline void check(string infoText)
  {
    #ifndef NDEBUG
    int rank;
    MPI_Comm comm = grid->getCartComm();
    MPI_Comm_rank(comm,&rank);
    MPI_Barrier(comm);

    #pragma omp parallel for schedule(static)
    for(int i=0; i<vInfo.size(); i++)
    {
      BlockInfo info = vInfo[i];
      FluidBlock& b = *(FluidBlock*)info.ptrBlock;

      for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
      for(int iy=0; iy<FluidBlock::sizeY; ++iy)
      for(int ix=0; ix<FluidBlock::sizeX; ++ix)
        if (std::isnan(b(ix,iy,iz).u) || std::isnan(b(ix,iy,iz).v) ||
            std::isnan(b(ix,iy,iz).w) || std::isnan(b(ix,iy,iz).p) )
          cout << infoText.c_str() << endl;
    }
    MPI_Barrier(comm);
    #endif
  }

  template <typename Kernel, bool skipBorder = false>
  void compute(const vector<Kernel*>& kernels)
  {
    SynchronizerMPI& Synch = grid->sync(*(kernels[0]));
    const int NX = grid->getBlocksPerDimension(0);
    const int NY = grid->getBlocksPerDimension(1);
    const int NZ = grid->getBlocksPerDimension(2);
    const int nthreads = omp_get_max_threads();
    LabMPI * labs = new LabMPI[nthreads];
    for(int i = 0; i < nthreads; ++i)
      labs[i].prepare(*grid, Synch);

    int rank;
    MPI_Comm_rank(grid->getCartComm(), &rank);
    MPI_Barrier(grid->getCartComm());
    vector<BlockInfo> avail0 = Synch.avail_inner();
    const int Ninner = avail0.size();

    #pragma omp parallel
    {
      int tid = omp_get_thread_num();
      Kernel& kernel=*(kernels[tid]);
      LabMPI& lab = labs[tid];

      #pragma omp for schedule(dynamic,1)
      for(int i=0; i<Ninner; i++) {
        BlockInfo I = avail0[i];
        if(skipBorder&& bSkip(NX, NY, NZ, I.index[0], I.index[1], I.index[2])){
          //printf("rank %d skipping inner %d %d %d %d %d %d\n", rank, I.index[0], I.index[1], I.index[2], NX, NY, NZ);
          continue;
        }

        FluidBlock& b = *(FluidBlock*)I.ptrBlock;
        lab.load(I, 0);
        kernel(lab, I, b);
      }
    }

    vector<BlockInfo> avail1 = Synch.avail_halo();
    const int Nhalo = avail1.size();

    #pragma omp parallel
    {
      int tid = omp_get_thread_num();
      Kernel& kernel=*(kernels[tid]);
      LabMPI& lab = labs[tid];

      #pragma omp for schedule(dynamic,1)
      for(int i=0; i<Nhalo; i++) {
        BlockInfo I = avail1[i];
        if(skipBorder&& bSkip(NX, NY, NZ, I.index[0], I.index[1], I.index[2])){
          //printf("rank %d skipping inner %d %d %d %d %d %d\n", rank, I.index[0], I.index[1], I.index[2], NX, NY, NZ);
          continue;
        }

        FluidBlock& b = *(FluidBlock*)I.ptrBlock;
        lab.load(I, 0);
        kernel(lab, I, b);
      }
    }

    if(labs!=NULL) {
      delete [] labs;
      labs=NULL;
    }

    MPI_Barrier(grid->getCartComm());
  }

  template <typename Kernel>
  void compute(const Kernel& kernel)
  {
    SynchronizerMPI& Synch = grid->sync(kernel);
    const int NX = grid->getBlocksPerDimension(0);
    const int NY = grid->getBlocksPerDimension(1);
    const int NZ = grid->getBlocksPerDimension(2);
    const int nthreads = omp_get_max_threads();
    LabMPI * labs = new LabMPI[nthreads];
    for(int i = 0; i < nthreads; ++i)
      labs[i].prepare(*grid, Synch);

    int rank;
    MPI_Comm_rank(grid->getCartComm(), &rank);
    MPI_Barrier(grid->getCartComm());
    vector<BlockInfo> avail0 = Synch.avail_inner();
    const int Ninner = avail0.size();

    #pragma omp parallel
    {
      int tid = omp_get_thread_num();
      LabMPI& lab = labs[tid];

      #pragma omp for schedule(dynamic,1)
      for(int i=0; i<Ninner; i++) {
        BlockInfo I = avail0[i];
        FluidBlock& b = *(FluidBlock*)I.ptrBlock;
        lab.load(I, 0);
        kernel(lab, I, b);
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
        BlockInfo I = avail1[i];
        FluidBlock& b = *(FluidBlock*)I.ptrBlock;
        lab.load(I, 0);
        kernel(lab, I, b);
      }
    }

    if(labs!=NULL) {
      delete [] labs;
      labs=NULL;
    }

    MPI_Barrier(grid->getCartComm());
  }

  static inline bool bSkip(const int nx, const int ny, const int nz, const int ix, const int iy, const int iz)
  {
    if(nx-1 == ix || 0 == ix) return true;
    if(ny-1 == iy || 0 == iy) return true;
    #ifndef BC_PERIODICZ
    if(nz-1 == iz || 0 == iz) return true;
    #endif
    return false;
  }

public:
  GenericCoordinator(FluidGridMPI * g) : grid(g)
  {
    vInfo = grid->getBlocksInfo();
  }
  virtual ~GenericCoordinator() {}
  virtual void operator()(const double dt) = 0;

  virtual string getName() = 0;
};

#endif
