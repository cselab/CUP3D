#pragma once

#include "Cubism/AMR_MeshAdaptationMPI.h"

using namespace cubism;

template <typename TGrid, typename TLab>
class MeshAdaptation_CUP : public MeshAdaptationMPI<TGrid,TLab>
{
 public:
   double Rtol_chi;
   double Ctol_chi;
   typedef typename TGrid::Block BlockType;
   typedef typename TGrid::BlockType::ElementType ElementType;
   typedef typename TGrid::BlockType Block;

   MeshAdaptation_CUP(TGrid &grid, double Rtol, double Ctol): MeshAdaptationMPI<TGrid,TLab>(grid,Rtol,Ctol)
   {
    Rtol_chi = 0.01;
    Ctol_chi = 0.001;
   }

   virtual void RefineBlocks(BlockType *B[8], BlockInfo parent) override
   {
      MeshAdaptationMPI<TGrid,TLab>::RefineBlocks(B,parent);

      //clipping for chi field
      for (int K = 0; K < 2; K++)
      for (int J = 0; J < 2; J++)
      for (int I = 0; I < 2; I++)
      {
        BlockType &b = *B[K * 4 + J * 2 + I];
        for (int k = 0; k < BlockType::sizeZ; k ++)
        for (int j = 0; j < BlockType::sizeY; j ++)
        for (int i = 0; i < BlockType::sizeX; i ++)
        {
          b(i  ,j  ,k  ).chi = std::max( b(i  ,j  ,k  ).chi, 0.0);
          b(i  ,j  ,k  ).chi = std::min( b(i  ,j  ,k  ).chi, 1.0);
        }
      }
   }

   virtual State TagLoadedBlock(TLab &Lab_, int level)
   {
      static const int nx = BlockType::sizeX;
      static const int ny = BlockType::sizeY;
      static const int nz = BlockType::sizeZ;
    
      double Linf   = 0.0; //chi
      double Linf_2 = 0.0; //vorticity
      for (int k = 0; k < nz; k++)
      for (int j = 0; j < ny; j++)
      for (int i = 0; i < nx; i++)
      {
        double s0 = std::fabs( Lab_(i, j, k).magnitude() );
        Linf = max(Linf,s0);
        const Real inv2h = .5;
        const ElementType &LW=Lab_(i-1,j  ,k  ), &LE=Lab_(i+1,j  ,k  );
        const ElementType &LS=Lab_(i  ,j-1,k  ), &LN=Lab_(i  ,j+1,k  );
        const ElementType &LF=Lab_(i  ,j  ,k-1), &LB=Lab_(i  ,j  ,k+1);
        double omega_x = inv2h * ( (LN.w-LS.w) - (LB.v-LF.v) );
        double omega_y = inv2h * ( (LB.u-LF.u) - (LE.w-LW.w) );
        double omega_z = inv2h * ( (LE.v-LW.v) - (LN.u-LS.u) );
        double omega = pow(omega_x*omega_x + omega_y*omega_y + omega_z*omega_z,0.5);
        if (i!=0 && j!=0 && k!=0 && i != nx-1 && j != ny-1 && k != nz-1)
          Linf_2 = max(Linf_2,omega);
      }
      //Linf   *= 1.0/(level+1);
      Linf_2 *= 1.0/(level+1);

      if (Linf > Rtol_chi || Linf_2 > MeshAdaptationMPI<TGrid,TLab>::tolerance_for_refinement ) return Refine;
      if (Linf < Ctol_chi && Linf_2 < MeshAdaptationMPI<TGrid,TLab>::tolerance_for_compression) return Compress;    
      return Leave;
   }
};