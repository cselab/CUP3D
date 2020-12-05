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
    Rtol_chi = 0.2;
    Ctol_chi = 0.000001;
   }

   virtual void RefineBlocks(BlockType *B[8], BlockInfo parent) override
   {
      static const int nx = BlockType::sizeX;
      static const int ny = BlockType::sizeY;
      static const int nz = BlockType::sizeZ;

#if 0 //derivatives
      MeshAdaptationMPI<TGrid,TLab>::RefineBlocks(B,parent);

#else //WENO3
      int tid = omp_get_thread_num();
      int offsetX[2] = {0, nx / 2};
      int offsetY[2] = {0, ny / 2};
      int offsetZ[2] = {0, nz / 2};
      TLab &Lab_ = MeshAdaptationMPI<TGrid,TLab>::labs[tid];

      for (int K = 0; K < 2; K++)
      for (int J = 0; J < 2; J++)
      for (int I = 0; I < 2; I++)
      {
        BlockType &b = *B[K * 4 + J * 2 + I];
        for (int k = 0; k < nz; k += 2)
        for (int j = 0; j < ny; j += 2)
        for (int i = 0; i < nx; i += 2)
        {
           const int Nweno = WENOWAVELET;
           ElementType El[Nweno][Nweno][Nweno];
           for (int i0 = -Nweno / 2; i0 <= Nweno / 2; i0++)
           for (int i1 = -Nweno / 2; i1 <= Nweno / 2; i1++)
           for (int i2 = -Nweno / 2; i2 <= Nweno / 2; i2++)
               El[i0 + Nweno / 2][i1 + Nweno / 2][i2 + Nweno / 2] = Lab_(i / 2 + offsetX[I] + i0,
                                                                         j / 2 + offsetY[J] + i1,
                                                                         k / 2 + offsetZ[K] + i2);

           ElementType Lines[Nweno][Nweno][2];
           ElementType Planes[Nweno][4];
           ElementType Ref[8];

           for (int i2 = -Nweno / 2; i2 <= Nweno / 2; i2++)
           for (int i1 = -Nweno / 2; i1 <= Nweno / 2; i1++)
                Kernel_1D(El[0][i1 + Nweno / 2][i2 + Nweno / 2],
                          El[1][i1 + Nweno / 2][i2 + Nweno / 2],
                          El[2][i1 + Nweno / 2][i2 + Nweno / 2],
                          Lines[i1 + Nweno / 2][i2 + Nweno / 2][0],
                          Lines[i1 + Nweno / 2][i2 + Nweno / 2][1]);
           for (int i2 = -Nweno / 2; i2 <= Nweno / 2; i2++)
           {
                Kernel_1D(Lines[0][i2 + Nweno / 2][0], Lines[1][i2 + Nweno / 2][0],
                          Lines[2][i2 + Nweno / 2][0], Planes[i2 + Nweno / 2][0],
                          Planes[i2 + Nweno / 2][1]);
                Kernel_1D(Lines[0][i2 + Nweno / 2][1], Lines[1][i2 + Nweno / 2][1],
                          Lines[2][i2 + Nweno / 2][1], Planes[i2 + Nweno / 2][2],
                          Planes[i2 + Nweno / 2][3]);
           }
           Kernel_1D(Planes[0][0], Planes[1][0], Planes[2][0], Ref[0], Ref[1]);
           Kernel_1D(Planes[0][1], Planes[1][1], Planes[2][1], Ref[2], Ref[3]);
           Kernel_1D(Planes[0][2], Planes[1][2], Planes[2][2], Ref[4], Ref[5]);
           Kernel_1D(Planes[0][3], Planes[1][3], Planes[2][3], Ref[6], Ref[7]);

           b(i, j, k)             = Ref[0];
           b(i, j, k + 1)         = Ref[1];
           b(i, j + 1, k)         = Ref[2];
           b(i, j + 1, k + 1)     = Ref[3];
           b(i + 1, j, k)         = Ref[4];
           b(i + 1, j, k + 1)     = Ref[5];
           b(i + 1, j + 1, k)     = Ref[6];
           b(i + 1, j + 1, k + 1) = Ref[7];
        }
      }
#endif
      //clipping for chi field
      for (int K = 0; K < 2; K++)
      for (int J = 0; J < 2; J++)
      for (int I = 0; I < 2; I++)
      {
        BlockType &b = *B[K * 4 + J * 2 + I];
        for (int iz = 0; iz < nz; iz ++)
        for (int iy = 0; iy < ny; iy ++)
        for (int ix = 0; ix < nx; ix ++)
        {
          b(ix  ,iy  ,iz  ).chi = std::max( b(ix  ,iy  ,iz  ).chi, 0.0);
          b(ix  ,iy  ,iz  ).chi = std::min( b(ix  ,iy  ,iz  ).chi, 1.0);
        }
      }
   }

   void Kernel_1D(ElementType E0, ElementType E1, ElementType E2, ElementType &left, ElementType &right)
   {
      for (int i = 0 ; i < ElementType::DIM ; i++)
        WENOWavelets3(E0.member(i), E1.member(i), E2.member(i), left.member(i), right.member(i));
   }
   void WENOWavelets3(double cm, double c, double cp, double &left, double &right)
   {
      double b1  = (c - cm) * (c - cm);
      double b2  = (c - cp) * (c - cp);
      double w1  = (1e-6 + b2) * (1e-6 + b2); // yes, 2 goes to 1 and 1 goes to 2
      double w2  = (1e-6 + b1) * (1e-6 + b1);
      double aux = 1.0 / (w1 + w2);
      w1 *= aux;
      w2 *= aux;
      double g1, g2;
      g1    = 0.75 * c + 0.25 * cm;
      g2    = 1.25 * c - 0.25 * cp;
      left  = g1 * w1 + g2 * w2;
      g1    = 1.25 * c - 0.25 * cm;
      g2    = 0.75 * c + 0.25 * cp;
      right = g1 * w1 + g2 * w2;
   }

   virtual State TagLoadedBlock(TLab &Lab_, int level, BlockInfo & info)
   {
      //This assumes zero Neumann BCs for velocity
      static const int nx = BlockType::sizeX;
      static const int ny = BlockType::sizeY;
      static const int nz = BlockType::sizeZ;
    
      double Linf   = 0.0; //chi
      double Linf_2 = 0.0; //vorticity

      for (int k = 0; k < nz; k++)
      for (int j = 0; j < ny; j++)
      for (int i = 0; i < nx; i++)
      {
        const Real inv2h = .5/info.h;
        const ElementType &LW=Lab_(i-1,j  ,k  ), &LE=Lab_(i+1,j  ,k  );
        const ElementType &LS=Lab_(i  ,j-1,k  ), &LN=Lab_(i  ,j+1,k  );
        const ElementType &LF=Lab_(i  ,j  ,k-1), &LB=Lab_(i  ,j  ,k+1);
        double omega_x = inv2h * ( (LN.w-LS.w) - (LB.v-LF.v) );
        double omega_y = inv2h * ( (LB.u-LF.u) - (LE.w-LW.w) );
        double omega_z = inv2h * ( (LE.v-LW.v) - (LN.u-LS.u) );
        double omega = pow(omega_x*omega_x + omega_y*omega_y + omega_z*omega_z,0.5);
        Linf_2 = max(Linf_2,omega);
      }

      for (int k = 0 - 1; k < nz + 1; k++)
      for (int j = 0 - 1; j < ny + 1; j++)
      for (int i = 0 - 1; i < nx + 1; i++)
      {
        double s0 = Lab_(i,j,k).chi;
        if (s0 > Linf) Linf = s0;
      }

      //Linf_2 *= 1.0/(level+1);

      if (Linf > Rtol_chi || Linf_2 > MeshAdaptationMPI<TGrid,TLab>::tolerance_for_refinement ) return Refine;
      if (Linf < Ctol_chi && Linf_2 < MeshAdaptationMPI<TGrid,TLab>::tolerance_for_compression) return Compress;    
      return Leave;
   }
};
