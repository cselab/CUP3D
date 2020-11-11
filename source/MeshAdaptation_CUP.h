#pragma once

#include "Cubism/AMR_MeshAdaptationMPI.h"

using namespace cubism;

template <typename TGrid, typename TLab>
class MeshAdaptation_CUP : public MeshAdaptationMPI<TGrid,TLab>
{
 public:
   double Rtol_chi;
   double Ctol_chi;

   MeshAdaptation_CUP(TGrid &grid, double Rtol, double Ctol): MeshAdaptationMPI<TGrid,TLab>(grid,Rtol,Ctol)
   {
    Rtol_chi = 0.01;
    Ctol_chi = 0.001;

   }

   virtual void RefineBlocks(BlockType *B[8], BlockInfo parent) override
   {
      int tid = omp_get_thread_num();
      const int nx = BlockType::sizeX;
      const int ny = BlockType::sizeY;

      int offsetX[2] = {0, nx / 2};
      int offsetY[2] = {0, ny / 2};

      TLab &Lab = labs[tid];

      const int nz = BlockType::sizeZ;
      int offsetZ[2] = {0, nz / 2};

      for (int K = 0; K < 2; K++)
      for (int J = 0; J < 2; J++)
      for (int I = 0; I < 2; I++)
      {
        BlockType &b = *B[K * 4 + J * 2 + I];
        b.clear();

        for (int k = 0; k < nz; k += 2)
        for (int j = 0; j < ny; j += 2)
        for (int i = 0; i < nx; i += 2)
        {
          ElementType dudx = 0.5*( Lab(i/2+offsetX[I]+1,j/2+offsetY[J]  ,k/2+offsetZ[K]  )-Lab(i/2+offsetX[I]-1,j/2+offsetY[J]  ,k/2+offsetZ[K]  ));
          ElementType dudy = 0.5*( Lab(i/2+offsetX[I]  ,j/2+offsetY[J]+1,k/2+offsetZ[K]  )-Lab(i/2+offsetX[I]  ,j/2+offsetY[J]-1,k/2+offsetZ[K]  ));
          ElementType dudz = 0.5*( Lab(i/2+offsetX[I]  ,j/2+offsetY[J]  ,k/2+offsetZ[K]+1)-Lab(i/2+offsetX[I]  ,j/2+offsetY[J]  ,k/2+offsetZ[K]-1));

          b(i  ,j  ,k  ) = Lab( i   /2+offsetX[I], j   /2+offsetY[J]  ,k    /2+offsetZ[K] )+ (2*( i   %2)-1)*0.25*dudx + (2*( j   %2)-1)*0.25*dudy + (2*(k    %2)-1)*0.25*dudz; 
          b(i+1,j  ,k  ) = Lab((i+1)/2+offsetX[I], j   /2+offsetY[J]  ,k    /2+offsetZ[K] )+ (2*((i+1)%2)-1)*0.25*dudx + (2*( j   %2)-1)*0.25*dudy + (2*(k    %2)-1)*0.25*dudz; 
          b(i  ,j+1,k  ) = Lab( i   /2+offsetX[I],(j+1)/2+offsetY[J]  ,k    /2+offsetZ[K] )+ (2*( i   %2)-1)*0.25*dudx + (2*((j+1)%2)-1)*0.25*dudy + (2*(k    %2)-1)*0.25*dudz; 
          b(i+1,j+1,k  ) = Lab((i+1)/2+offsetX[I],(j+1)/2+offsetY[J]  ,k    /2+offsetZ[K] )+ (2*((i+1)%2)-1)*0.25*dudx + (2*((j+1)%2)-1)*0.25*dudy + (2*(k    %2)-1)*0.25*dudz; 
          b(i  ,j  ,k+1) = Lab( i   /2+offsetX[I], j   /2+offsetY[J]  ,(k+1)/2+offsetZ[K] )+ (2*( i   %2)-1)*0.25*dudx + (2*( j   %2)-1)*0.25*dudy + (2*((k+1)%2)-1)*0.25*dudz; 
          b(i+1,j  ,k+1) = Lab((i+1)/2+offsetX[I], j   /2+offsetY[J]  ,(k+1)/2+offsetZ[K] )+ (2*((i+1)%2)-1)*0.25*dudx + (2*( j   %2)-1)*0.25*dudy + (2*((k+1)%2)-1)*0.25*dudz; 
          b(i  ,j+1,k+1) = Lab( i   /2+offsetX[I],(j+1)/2+offsetY[J]  ,(k+1)/2+offsetZ[K] )+ (2*( i   %2)-1)*0.25*dudx + (2*((j+1)%2)-1)*0.25*dudy + (2*((k+1)%2)-1)*0.25*dudz; 
          b(i+1,j+1,k+1) = Lab((i+1)/2+offsetX[I],(j+1)/2+offsetY[J]  ,(k+1)/2+offsetZ[K] )+ (2*((i+1)%2)-1)*0.25*dudx + (2*((j+1)%2)-1)*0.25*dudy + (2*((k+1)%2)-1)*0.25*dudz;

          //clipping for chi field
          b(i  ,j  ,k  ).chi = std::min( b(i  ,j  ,k  ).chi, 0.0);
          b(i+1,j  ,k  ).chi = std::min( b(i+1,j  ,k  ).chi, 0.0);
          b(i  ,j+1,k  ).chi = std::min( b(i  ,j+1,k  ).chi, 0.0);
          b(i+1,j+1,k  ).chi = std::min( b(i+1,j+1,k  ).chi, 0.0);
          b(i  ,j  ,k+1).chi = std::min( b(i  ,j  ,k+1).chi, 0.0);
          b(i+1,j  ,k+1).chi = std::min( b(i+1,j  ,k+1).chi, 0.0);
          b(i  ,j+1,k+1).chi = std::min( b(i  ,j+1,k+1).chi, 0.0);
          b(i+1,j+1,k+1).chi = std::min( b(i+1,j+1,k+1).chi, 0.0);

          b(i  ,j  ,k  ).chi = std::max( b(i  ,j  ,k  ).chi, 1.0);
          b(i+1,j  ,k  ).chi = std::max( b(i+1,j  ,k  ).chi, 1.0);
          b(i  ,j+1,k  ).chi = std::max( b(i  ,j+1,k  ).chi, 1.0);
          b(i+1,j+1,k  ).chi = std::max( b(i+1,j+1,k  ).chi, 1.0);
          b(i  ,j  ,k+1).chi = std::max( b(i  ,j  ,k+1).chi, 1.0);
          b(i+1,j  ,k+1).chi = std::max( b(i+1,j  ,k+1).chi, 1.0);
          b(i  ,j+1,k+1).chi = std::max( b(i  ,j+1,k+1).chi, 1.0);
          b(i+1,j+1,k+1).chi = std::max( b(i+1,j+1,k+1).chi, 1.0);
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
      }
      Linf   *= 1.0/(level+1);
      Linf_2 *= 1.0/(level+1);
      if (Linf > Rtol_chi || Linf_2 > MeshAdaptationMPI<TGrid,TLab>::tolerance_for_refinement ) return Refine;
      if (Linf < Ctol_chi && Linf_2 < MeshAdaptationMPI<TGrid,TLab>::tolerance_for_compression) return Compress;    
      return Leave;
   }
};