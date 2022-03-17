//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//

#pragma once
#include "../SimulationData.h"
#include <Cubism/FluxCorrectionMPI.h>

CubismUP_3D_NAMESPACE_BEGIN

class Operator
{
  protected:
   SimulationData & sim;

   template <typename Kernel, typename TGrid, typename Lab, typename TGrid_corr = TGrid>
   void compute(Kernel & kernel, TGrid * g, TGrid_corr * g_corr = nullptr)
   {
     if (g_corr != nullptr) g_corr->Corrector.prepare(*g_corr);

     cubism::SynchronizerMPI_AMR<typename TGrid::Real,TGrid>& Synch = *(g->sync(kernel));
     std::vector<cubism::BlockInfo*> *inner = &Synch.avail_inner();
     std::vector<cubism::BlockInfo*> *halo;
     #pragma omp parallel
     {
         Lab lab;
         lab.prepare(*g, Synch);

         #pragma omp for nowait
         for (const cubism::BlockInfo *I : *inner)
         {
           lab.load(*I, 0);
           kernel(lab, *I);
         }

         #pragma omp master
         halo = &Synch.avail_halo();
         #pragma omp barrier

         #pragma omp for nowait
         for (const cubism::BlockInfo *I : *halo)
         {
           lab.load(*I, 0);
           kernel(lab, *I);
         }
     }

     if (g_corr != nullptr) g_corr->Corrector.FillBlockCases();
   }

   template <typename Kernel, typename TGrid, typename Lab, typename TGrid_corr = TGrid>
   void compute(const Kernel & kernel, TGrid * g, TGrid_corr * g_corr = nullptr)
   {
     if (g_corr != nullptr) g_corr->Corrector.prepare(*g_corr);

     cubism::SynchronizerMPI_AMR<typename TGrid::Real,TGrid>& Synch = *(g->sync(kernel));
     std::vector<cubism::BlockInfo*> *inner = &Synch.avail_inner();
     std::vector<cubism::BlockInfo*> *halo;
     #pragma omp parallel
     {
         Lab lab;
         lab.prepare(*g, Synch);

         #pragma omp for nowait
         for (const cubism::BlockInfo *I : *inner)
         {
           lab.load(*I, 0);
           kernel(lab, *I);
         }

         #pragma omp master
         halo = &Synch.avail_halo();
         #pragma omp barrier

         #pragma omp for nowait
         for (const cubism::BlockInfo *I : *halo)
         {
           lab.load(*I, 0);
           kernel(lab, *I);
         }
     }

     if (g_corr != nullptr) g_corr->Corrector.FillBlockCases();
   }

  public:
   Operator(SimulationData & s) : sim(s) {  }
   virtual ~Operator() = default;
   virtual void operator()(const Real dt) = 0;
   virtual std::string getName() = 0;
};
CubismUP_3D_NAMESPACE_END
