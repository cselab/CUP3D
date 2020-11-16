//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch) and Christian Conti.
//

#ifndef CubismUP_3D_Operator_h
#define CubismUP_3D_Operator_h

#include "../SimulationData.h"

CubismUP_3D_NAMESPACE_BEGIN

class Operator
{
 protected:
  SimulationData & sim;
  FluidGridMPI * const grid = sim.grid;

  inline void check(const std::string &infoText)
  {
    #ifndef NDEBUG
    std::vector<cubism::BlockInfo>& vInfo = sim.vInfo();
    MPI_Comm comm = grid->getCartComm();
    MPI_Barrier(comm);

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < (int)vInfo.size(); ++i)
    {
      const cubism::BlockInfo info = vInfo[i];
      const FluidBlock& b = *(FluidBlock*)info.ptrBlock;
      for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
      for(int iy=0; iy<FluidBlock::sizeY; ++iy)
      for(int ix=0; ix<FluidBlock::sizeX; ++ix)
        if (std::isnan(b(ix,iy,iz).u   ) || std::isnan(b(ix,iy,iz).v   ) ||
            std::isnan(b(ix,iy,iz).w   ) || std::isnan(b(ix,iy,iz).p   ) ||
            std::isnan(b(ix,iy,iz).tmpU) || std::isnan(b(ix,iy,iz).tmpV) ||
            std::isnan(b(ix,iy,iz).tmpW) || std::isnan(b(ix,iy,iz).chi) )
        {
          fflush(stderr);
          std::cout << "*************" << std::endl;
          std::cout << "BLOCK (" << info.index[0] << "," << info.index[1] << "," << info.index[2] << ")" << std::endl;
          std::cout << "level = " << info.level << std::endl;
          std::cout << "Z = " << info.Z << std::endl;
          std::cout << "ix=" << ix << std::endl;
          std::cout << "iy=" << iy << std::endl;
          std::cout << "iz=" << iz << std::endl;
          std::cout << "u =" << b(ix,iy,iz).u << std::endl;
          std::cout << "v =" << b(ix,iy,iz).v << std::endl;
          std::cout << "w =" << b(ix,iy,iz).w << std::endl;
          std::cout << "p =" << b(ix,iy,iz).p << std::endl;
          std::cout << "tmpU =" << b(ix,iy,iz).tmpU << std::endl;
          std::cout << "tmpV =" << b(ix,iy,iz).tmpV << std::endl;
          std::cout << "tmpW =" << b(ix,iy,iz).tmpW << std::endl;
          std::cout << "chi =" << b(ix,iy,iz).chi << std::endl;
          printf("GenericCoordinator::check isnan %s\n", infoText.c_str());
          std::cout << "*************" << std::endl;
          fflush(stdout); MPI_Abort(comm, 1);
        }
    }
    MPI_Barrier(comm);
    #endif
  }

  template <typename Kernel>
  void compute(const std::vector<Kernel*>& kernels)
  {
    cubism::SynchronizerMPI_AMR<Real,FluidGridMPI>& Synch = *grid->sync(*(kernels[0]));
    const int nthreads = omp_get_max_threads();
    LabMPI * labs = new LabMPI[nthreads];
    #pragma omp parallel for schedule(static, 1)
    for(int i = 0; i < nthreads; ++i) {
      labs[i].setBC(sim.BCx_flag, sim.BCy_flag, sim.BCz_flag);
      labs[i].prepare(* sim.grid, Synch);
    }

    int rank;
    MPI_Comm_rank(grid->getCartComm(), &rank);
    MPI_Barrier(grid->getCartComm());
    std::vector<cubism::BlockInfo*> avail0 = Synch.avail_inner();
    const int Ninner = avail0.size();

    #pragma omp parallel
    {
      int tid = omp_get_thread_num();
      Kernel& kernel = * (kernels[tid]); LabMPI& lab = labs[tid];

      #pragma omp for schedule(static)
      for(int i=0; i<Ninner; i++) {
        const cubism::BlockInfo& I = *avail0[i];
        FluidBlock& b = *(FluidBlock*)I.ptrBlock;
        lab.load(I, 0);
        kernel(lab, I, b);
      }
    }

    if(sim.nprocs>1)
    {
      std::vector<cubism::BlockInfo*> avail1 = Synch.avail_halo();
      const int Nhalo = avail1.size();

      #pragma omp parallel
      {
        int tid = omp_get_thread_num();
        Kernel& kernel = * (kernels[tid]); LabMPI& lab = labs[tid];

        #pragma omp for schedule(static)
        for(int i=0; i<Nhalo; i++) {
          const cubism::BlockInfo& I = *avail1[i];
          FluidBlock& b = *(FluidBlock*)I.ptrBlock;
          lab.load(I, 0);
          kernel(lab, I, b);
        }
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
    cubism::SynchronizerMPI_AMR<Real,FluidGridMPI>& Synch = *grid->sync(kernel);

    const int nthreads = omp_get_max_threads();
    LabMPI * labs = new LabMPI[nthreads];
    #pragma omp parallel for schedule(static, 1)
    for(int i = 0; i < nthreads; ++i) {
      labs[i].setBC(sim.BCx_flag, sim.BCy_flag, sim.BCz_flag);
      labs[i].prepare(* sim.grid, Synch);
    }

    int rank;
    MPI_Comm_rank(grid->getCartComm(), &rank);
    MPI_Barrier(grid->getCartComm());
    std::vector<cubism::BlockInfo*> avail0 = Synch.avail_inner();
    const int Ninner = avail0.size();

    #pragma omp parallel
    {
      int tid = omp_get_thread_num();
      LabMPI& lab = labs[tid];

      #pragma omp for schedule(static)
      for(int i=0; i<Ninner; i++) {
        const cubism::BlockInfo I = *avail0[i];
        FluidBlock& b = *(FluidBlock*)I.ptrBlock;
        lab.load(I, 0);
        kernel(lab, I, b);
      }
    }

    if(sim.nprocs>1)
    {
      std::vector<cubism::BlockInfo*> avail1 = Synch.avail_halo();
      const int Nhalo = avail1.size();

      #pragma omp parallel
      {
        int tid = omp_get_thread_num();
        LabMPI& lab = labs[tid];

        #pragma omp for schedule(static)
        for(int i=0; i<Nhalo; i++) {
          const cubism::BlockInfo I = *avail1[i];
          FluidBlock& b = *(FluidBlock*)I.ptrBlock;
          lab.load(I, 0);
          kernel(lab, I, b);
        }
      }
    }

    if(labs != nullptr) {
      delete [] labs;
      labs = nullptr;
    }

    MPI_Barrier(grid->getCartComm());
  }

public:
  Operator(SimulationData & s) : sim(s) {  }
  virtual ~Operator() = default;
  virtual void operator()(const double dt) = 0;
  virtual std::string getName() = 0;
};

CubismUP_3D_NAMESPACE_END
#endif
