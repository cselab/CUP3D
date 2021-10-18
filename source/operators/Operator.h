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
#include <Cubism/FluxCorrectionMPI.h>

CubismUP_3D_NAMESPACE_BEGIN

class Operator
{
 protected:
  SimulationData & sim;
  FluidGridMPI * const grid = sim.grid;
  FluidGridMPIPoisson * const gridPoisson = sim.gridPoisson;

  inline void check(const std::string &infoText)
  {
    #ifndef NDEBUG
    std::vector<cubism::BlockInfo>& vInfo = sim.vInfo();
    MPI_Comm comm = grid->getCartComm();
    MPI_Barrier(comm);

    //#pragma omp parallel for schedule(static)
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

    if (sim.step > sim.step_2nd_start)
    for (int i = 0; i < (int)vInfo.size(); ++i)
    {
      const cubism::BlockInfo info = vInfo[i];
      const FluidBlock& b = *(FluidBlock*)info.ptrBlock;
      for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
      for(int iy=0; iy<FluidBlock::sizeY; ++iy)
      for(int ix=0; ix<FluidBlock::sizeX; ++ix)
        if (std::isnan(b.dataOld[iz][iy][ix][0]) ||
            std::isnan(b.dataOld[iz][iy][ix][1]) ||
            std::isnan(b.dataOld[iz][iy][ix][2]) )
        {
          fflush(stderr);
          std::cout << "*************" << std::endl;
          std::cout << "BLOCK (" << info.index[0] << "," << info.index[1] << "," << info.index[2] << ")" << std::endl;
          std::cout << "level = " << info.level << std::endl;
          std::cout << "Z = " << info.Z << std::endl;
          std::cout << "ix=" << ix << std::endl;
          std::cout << "iy=" << iy << std::endl;
          std::cout << "iz=" << iz << std::endl;
          std::cout << "data u =" << b.dataOld[iz][iy][ix][0] << std::endl;
          std::cout << "data v =" << b.dataOld[iz][iy][ix][1] << std::endl;
          std::cout << "data w =" << b.dataOld[iz][iy][ix][2] << std::endl;
          std::cout << "data p =" << b.dataOld[iz][iy][ix][3] << std::endl;
          printf("GenericCoordinator::check data isnan %s\n", infoText.c_str());
          std::cout << "*************" << std::endl;
          fflush(stdout); MPI_Abort(comm, 1);
        }
    }
    MPI_Barrier(comm);
    #endif
  }

  /// Execute kernel for each block, where a different kernel instance is
  /// provided for each thread.
  template <typename Kernel>
  void compute(const std::vector<Kernel*>& kernels, const bool applyFluxCorrection = false, const bool FluxIntegration = false)
  {
    if (applyFluxCorrection) sim.grid->Corrector.prepare(*sim.grid);

    cubism::SynchronizerMPI_AMR<Real,FluidGridMPI>& Synch = *grid->sync(*(kernels[0]));
    const int nthreads = omp_get_max_threads();
    LabMPI * labs = new LabMPI[nthreads];
    #pragma omp parallel for schedule(static, 1)
    for(int i = 0; i < nthreads; ++i) {
      labs[i].setBC(sim.BCx_flag, sim.BCy_flag, sim.BCz_flag);
      labs[i].prepare(* sim.grid, Synch);
    }

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

    if (applyFluxCorrection) sim.grid->Corrector.FillBlockCases(FluxIntegration);
  }

  template <typename Kernel>
  void compute(const Kernel& kernel, const bool applyFluxCorrection = false, const bool FluxIntegration = false)
  {
    if (applyFluxCorrection) sim.grid->Corrector.prepare(*sim.grid);
    cubism::SynchronizerMPI_AMR<Real,FluidGridMPI>& Synch = *grid->sync(kernel);

    const int nthreads = omp_get_max_threads();
    LabMPI * labs = new LabMPI[nthreads];
    #pragma omp parallel for schedule(static, 1)
    for(int i = 0; i < nthreads; ++i) {
      labs[i].setBC(sim.BCx_flag, sim.BCy_flag, sim.BCz_flag);
      labs[i].prepare(* sim.grid, Synch);
    }

    std::vector<cubism::BlockInfo*> avail0 = Synch.avail_inner();
    const int Ninner = avail0.size();

    #pragma omp parallel
    {
      int tid = omp_get_thread_num();
      LabMPI& lab = labs[tid];

      #pragma omp for schedule(static)
      for(int i=0; i<Ninner; i++) {
        const cubism::BlockInfo &I = *avail0[i];
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
          const cubism::BlockInfo &I = *avail1[i];
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

    if (applyFluxCorrection) sim.grid->Corrector.FillBlockCases(FluxIntegration);
  }

  template <typename Kernel>
  void computePoisson(const Kernel& kernel, const bool applyFluxCorrection = false, const bool FluxIntegration = false)
  {
    if (applyFluxCorrection) sim.gridPoisson->Corrector.prepare(*sim.gridPoisson);

    cubism::SynchronizerMPI_AMR<Real,FluidGridMPIPoisson>& Synch = *gridPoisson->sync(kernel);

    const int nthreads = omp_get_max_threads();
    LabMPIPoisson * labs = new LabMPIPoisson[nthreads];
    #pragma omp parallel for schedule(static, 1)
    for(int i = 0; i < nthreads; ++i) {
      labs[i].setBC(sim.BCx_flag, sim.BCy_flag, sim.BCz_flag);
      labs[i].prepare(* sim.gridPoisson, Synch);
    }

    std::vector<cubism::BlockInfo*> avail0 = Synch.avail_inner();
    const int Ninner = avail0.size();


    #pragma omp parallel
    {
      int tid = omp_get_thread_num();
      LabMPIPoisson& lab = labs[tid];

      #pragma omp for
      for(int i=0; i<Ninner; i++) {
        const cubism::BlockInfo& I = *avail0[i];
        FluidBlockPoisson& b = *(FluidBlockPoisson*)I.ptrBlock;
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
        LabMPIPoisson& lab = labs[tid];

        #pragma omp for
        for(int i=0; i<Nhalo; i++) {
          const cubism::BlockInfo& I = *avail1[i];
          FluidBlockPoisson& b = *(FluidBlockPoisson*)I.ptrBlock;
          lab.load(I, 0);
          kernel(lab, I, b);
        }
      }
    }

    if(labs != nullptr) {
      delete [] labs;
      labs = nullptr;
    }
    if (applyFluxCorrection) sim.gridPoisson->Corrector.FillBlockCases(FluxIntegration);
  }

public:
  Operator(SimulationData & s) : sim(s) {  }
  virtual ~Operator() = default;
  virtual void operator()(const double dt) = 0;
  virtual std::string getName() = 0;
};

CubismUP_3D_NAMESPACE_END
#endif
