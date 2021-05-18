//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch) and Christian Conti.
//

#ifndef CubismUP_3D_ProcessOperators_h
#define CubismUP_3D_ProcessOperators_h

#include "../SimulationData.h"
#include "Operator.h"

CubismUP_3D_NAMESPACE_BEGIN

#ifndef CUP_SINGLE_PRECISION
#define MPIREAL MPI_DOUBLE
#else
#define MPIREAL MPI_FLOAT
#endif /* CUP_SINGLE_PRECISION */

inline Real findMaxUzeroMom(SimulationData& sim)
{
  const std::vector<cubism::BlockInfo>& myInfo = sim.vInfo();
  const Real uinf[3] = {sim.uinf[0], sim.uinf[1], sim.uinf[2]};
  Real mom[3] = {0, 0, 0};
  #pragma omp parallel for schedule(static) reduction(+ : mom[:3])
  for(size_t i=0; i<myInfo.size(); i++)
  {
    const cubism::BlockInfo& info = myInfo[i];
    const FluidBlock& b = *(const FluidBlock *)info.ptrBlock;

    for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
    for(int iy=0; iy<FluidBlock::sizeY; ++iy)
    for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
      Real h[3]; info.spacing(h, ix, iy, iz);
      const Real vol = h[0] * h[1] * h[2];
      mom[0] += vol * b(ix,iy,iz).u;
      mom[1] += vol * b(ix,iy,iz).v;
      mom[2] += vol * b(ix,iy,iz).w;
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, mom, 3, MPIREAL, MPI_SUM, sim.app_comm);
  const Real corrX = mom[0] / (sim.extent[0] * sim.extent[1] * sim.extent[2]);
  const Real corrY = mom[1] / (sim.extent[0] * sim.extent[1] * sim.extent[2]);
  const Real corrZ = mom[2] / (sim.extent[0] * sim.extent[1] * sim.extent[2]);
  if(sim.verbose)
    printf("Correction in relative momenta:[%e %e %e]\n",corrX,corrY,corrZ);

  Real maxU = 0;
  #pragma omp parallel for schedule(static) reduction(max : maxU)
  for(size_t i=0; i<myInfo.size(); i++)
  {
    const cubism::BlockInfo& info = myInfo[i];
    FluidBlock& b = *(FluidBlock *)info.ptrBlock;

    for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
    for(int iy=0; iy<FluidBlock::sizeY; ++iy)
    for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
      b(ix,iy,iz).u -= corrX; const Real u = std::fabs(b(ix,iy,iz).u +uinf[0]);
      b(ix,iy,iz).v -= corrY; const Real v = std::fabs(b(ix,iy,iz).v +uinf[1]);
      b(ix,iy,iz).w -= corrZ; const Real w = std::fabs(b(ix,iy,iz).w +uinf[2]);
      const Real maxUabsAdv = std::max({u, v, w});
      maxU = std::max(maxU, maxUabsAdv);
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, & maxU, 1, MPIREAL, MPI_MAX, sim.app_comm);
  assert(maxU >= 0);
  return maxU;
}

inline Real findMaxU(SimulationData& sim)
{
  const std::vector<cubism::BlockInfo>& myInfo = sim.vInfo();
  const Real uinf[3] = {sim.uinf[0], sim.uinf[1], sim.uinf[2]};

  Real maxU = 0;
  #pragma omp parallel for schedule(static) reduction(max : maxU)
  for(size_t i=0; i<myInfo.size(); i++)
  {
    const cubism::BlockInfo& info = myInfo[i];
    const FluidBlock& b = *(const FluidBlock *)info.ptrBlock;

    for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
    for(int iy=0; iy<FluidBlock::sizeY; ++iy)
    for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
      const Real advu = std::fabs(b(ix,iy,iz).u + uinf[0]);
      const Real advv = std::fabs(b(ix,iy,iz).v + uinf[1]);
      const Real advw = std::fabs(b(ix,iy,iz).w + uinf[2]);
      const Real maxUl = std::max({advu, advv, advw});
      maxU = std::max(maxU, maxUl);
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, & maxU, 1, MPIREAL, MPI_MAX, sim.app_comm);
  assert(maxU >= 0);
  return maxU;
}

#undef MPIREAL

class KernelVorticity
{
  public:
  KernelVorticity() = default;
  const std::array<int, 3> stencil_start = {-1,-1,-1}, stencil_end = {2, 2, 2};
  const cubism::StencilInfo stencil{-1,-1,-1, 2,2,2, false, {FE_U,FE_V,FE_W}};

  template <typename Lab, typename BlockType>
  void operator()(Lab & lab, const cubism::BlockInfo& info, BlockType& o) const
  {
    const Real inv2h = .5 * info.h_gridpoint * info.h_gridpoint;
    for (int iz=0; iz<FluidBlock::sizeZ; ++iz)
    for (int iy=0; iy<FluidBlock::sizeY; ++iy)
    for (int ix=0; ix<FluidBlock::sizeX; ++ix) {
      const FluidElement &LW=lab(ix-1,iy,iz), &LE=lab(ix+1,iy,iz);
      const FluidElement &LS=lab(ix,iy-1,iz), &LN=lab(ix,iy+1,iz);
      const FluidElement &LF=lab(ix,iy,iz-1), &LB=lab(ix,iy,iz+1);
      o(ix,iy,iz).tmpU = inv2h * ( (LN.w-LS.w) - (LB.v-LF.v) );
      o(ix,iy,iz).tmpV = inv2h * ( (LB.u-LF.u) - (LE.w-LW.w) );
      o(ix,iy,iz).tmpW = inv2h * ( (LE.v-LW.v) - (LN.u-LS.u) );
    }
    BlockCase<BlockType> * tempCase = (BlockCase<BlockType> *)(info.auxiliary);
    typename BlockType::ElementType * faceXm = nullptr;
    typename BlockType::ElementType * faceXp = nullptr;
    typename BlockType::ElementType * faceYm = nullptr;
    typename BlockType::ElementType * faceYp = nullptr;
    typename BlockType::ElementType * faceZp = nullptr;
    typename BlockType::ElementType * faceZm = nullptr;
    if (tempCase != nullptr)
    {
        faceXm = tempCase -> storedFace[0] ?  & tempCase -> m_pData[0][0] : nullptr;
        faceXp = tempCase -> storedFace[1] ?  & tempCase -> m_pData[1][0] : nullptr;
        faceYm = tempCase -> storedFace[2] ?  & tempCase -> m_pData[2][0] : nullptr;
        faceYp = tempCase -> storedFace[3] ?  & tempCase -> m_pData[3][0] : nullptr;
        faceZm = tempCase -> storedFace[4] ?  & tempCase -> m_pData[4][0] : nullptr;
        faceZp = tempCase -> storedFace[5] ?  & tempCase -> m_pData[5][0] : nullptr;
    }
    if (faceXm != nullptr)
    {
        int ix = 0;
        for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
        for(int iy=0; iy<FluidBlock::sizeY; ++iy)
        {
          const FluidElement &LW=lab(ix-1,iy,iz);
          const FluidElement &LC=lab(ix,iy,iz);
          faceXm[iy + FluidBlock::sizeY * iz].clear();
          faceXm[iy + FluidBlock::sizeY * iz].tmpV = -inv2h*( LW.w + LC.w );
          faceXm[iy + FluidBlock::sizeY * iz].tmpW = +inv2h*( LW.v + LC.v );
        }
    }
    if (faceXp != nullptr)
    {
       int ix = FluidBlock::sizeX-1;
       for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
       for(int iy=0; iy<FluidBlock::sizeY; ++iy)
       {
          const FluidElement &LE=lab(ix+1,iy,iz);
          const FluidElement &LC=lab(ix,iy,iz);
          faceXp[iy + FluidBlock::sizeY * iz].clear();
          faceXp[iy + FluidBlock::sizeY * iz].tmpV = +inv2h*( LE.w + LC.w );
          faceXp[iy + FluidBlock::sizeY * iz].tmpW = -inv2h*( LE.v + LC.v );
       }
    }
    if (faceYm != nullptr)
    {
       int iy = 0;
       for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
       for(int ix=0; ix<FluidBlock::sizeX; ++ix)
       {
          const FluidElement &LS=lab(ix,iy-1,iz);
          const FluidElement &LC=lab(ix,iy,iz);
          faceYm[ix + FluidBlock::sizeX * iz].clear();
          faceYm[ix + FluidBlock::sizeX * iz].tmpU = +inv2h*(LS.w+LC.w);
          faceYm[ix + FluidBlock::sizeX * iz].tmpW = -inv2h*(LS.u+LC.u);
       }
     }
     if (faceYp != nullptr)
     {
       int iy = FluidBlock::sizeY-1;
       for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
       for(int ix=0; ix<FluidBlock::sizeX; ++ix)
       {
          const FluidElement &LN=lab(ix,iy+1,iz);
          const FluidElement &LC=lab(ix,iy,iz);
          faceYp[ix + FluidBlock::sizeX * iz].clear();
          faceYp[ix + FluidBlock::sizeX * iz].tmpU = -inv2h*(LN.w+LC.w);
          faceYp[ix + FluidBlock::sizeX * iz].tmpW = +inv2h*(LN.u+LC.u);
       }
     }
     if (faceZm != nullptr)
     {
       int iz = 0;
       for(int iy=0; iy<FluidBlock::sizeY; ++iy)
       for(int ix=0; ix<FluidBlock::sizeX; ++ix)
       {
          const FluidElement &LF=lab(ix,iy,iz-1);
          const FluidElement &LC=lab(ix,iy,iz);
          faceZm[ix + FluidBlock::sizeX * iy].clear();
          faceZm[ix + FluidBlock::sizeX * iy].tmpU = -inv2h*(LF.v+LC.v);
          faceZm[ix + FluidBlock::sizeX * iy].tmpV = +inv2h*(LF.u+LC.u);
       }
     }
     if (faceZp != nullptr)
     {
       int iz = FluidBlock::sizeZ-1;
       for(int iy=0; iy<FluidBlock::sizeY; ++iy)
       for(int ix=0; ix<FluidBlock::sizeX; ++ix)
       {
          const FluidElement &LB=lab(ix,iy,iz+1);
          const FluidElement &LC=lab(ix,iy,iz);
          faceZp[ix + FluidBlock::sizeX * iy].clear();
          faceZp[ix + FluidBlock::sizeX * iy].tmpU = +inv2h*(LB.v+LC.v);
          faceZp[ix + FluidBlock::sizeX * iy].tmpV = -inv2h*(LB.u+LC.u);
       }
     }
  }
};

class ComputeVorticity : public Operator
{
  public:
  ComputeVorticity(SimulationData & s) : Operator(s) { }
  bool zeroX{false};
  bool zeroY{false};
  bool zeroZ{false};
  void operator()(const double dt)
  {
    sim.startProfiler("Vorticity Kernel");
    const KernelVorticity K;
    compute<KernelVorticity>(K,true);
    const std::vector<cubism::BlockInfo>& myInfo = sim.vInfo();
    #pragma omp parallel for
    for(size_t i=0; i<myInfo.size(); i++)
    {
      const cubism::BlockInfo& info = myInfo[i];
      FluidBlock& b = *( FluidBlock *)info.ptrBlock;
      const double fac = 1.0/(info.h*info.h*info.h);
      for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
      for(int iy=0; iy<FluidBlock::sizeY; ++iy)
      for(int ix=0; ix<FluidBlock::sizeX; ++ix)
      {
        b(ix,iy,iz).tmpU *=fac;
        b(ix,iy,iz).tmpV *=fac;
        b(ix,iy,iz).tmpW *=fac;
      }
    }

    if (zeroX || zeroY || zeroZ)
    {
      #pragma omp parallel for
      for(size_t i=0; i<myInfo.size(); i++)
      {
        const cubism::BlockInfo& info = myInfo[i];
        FluidBlock& b = *( FluidBlock *)info.ptrBlock;
        for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
        for(int iy=0; iy<FluidBlock::sizeY; ++iy)
        for(int ix=0; ix<FluidBlock::sizeX; ++ix)
        {
          if (zeroX) b(ix,iy,iz).tmpU = 0;
          if (zeroY) b(ix,iy,iz).tmpV = 0;
          if (zeroZ) b(ix,iy,iz).tmpW = 0;
        }
      }
    }
    sim.stopProfiler();
    check("Vorticity");
  }
  std::string getName() { return "Vorticity"; }
};

class KernelQcriterion
{
  public:
  KernelQcriterion() = default;
  const std::array<int, 3> stencil_start = {-1,-1,-1}, stencil_end = {2, 2, 2};
  const cubism::StencilInfo stencil{-1,-1,-1, 2,2,2, false, {FE_U,FE_V,FE_W}};

  template <typename Lab, typename BlockType>
  void operator()(Lab & lab, const cubism::BlockInfo& info, BlockType& o) const
  {
    const Real inv2h = .5 / info.h_gridpoint;
    for (int iz=0; iz<FluidBlock::sizeZ; ++iz)
    for (int iy=0; iy<FluidBlock::sizeY; ++iy)
    for (int ix=0; ix<FluidBlock::sizeX; ++ix) {
      const FluidElement &LW=lab(ix-1,iy,iz), &LE=lab(ix+1,iy,iz);
      const FluidElement &LS=lab(ix,iy-1,iz), &LN=lab(ix,iy+1,iz);
      const FluidElement &LF=lab(ix,iy,iz-1), &LB=lab(ix,iy,iz+1);
      const Real WX  = inv2h * ( (LN.w-LS.w) - (LB.v-LF.v) );
      const Real WY  = inv2h * ( (LB.u-LF.u) - (LE.w-LW.w) );
      const Real WZ  = inv2h * ( (LE.v-LW.v) - (LN.u-LS.u) );
      const Real D11 = inv2h * (LE.u-LW.u); // shear stresses
      const Real D22 = inv2h * (LN.v-LS.v); // shear stresses
      const Real D33 = inv2h * (LB.w-LF.w); // shear stresses
      const Real D12 = inv2h * (LN.u-LS.u + LE.v-LW.v); // shear stresses
      const Real D13 = inv2h * (LE.w-LW.w + LB.u-LF.u); // shear stresses
      const Real D23 = inv2h * (LB.v-LF.v + LN.w-LS.w); // shear stresses
      // trace( S S^t ) where S is the sym part of the vel gradient:
      const Real SS = D11*D11 +D22*D22 +D33*D33 +(D12*D12 +D13*D13 +D23*D23)/2;
      o(ix,iy,iz).p = ( (WX*WX + WY*WY + WZ*WZ)/2 - SS ) / 2;
    }
  }
};

class ComputeQcriterion : public Operator
{
  public:
  ComputeQcriterion(SimulationData & s) : Operator(s) { }
  void operator()(const double dt)
  {
    sim.startProfiler("Qcriterion Kernel");
    const KernelQcriterion K;
    compute<KernelQcriterion>(K);
    sim.stopProfiler();
    check("Qcriterion");
  }
  std::string getName() { return "Qcriterion"; }
};

class KernelDivergence
{
  public:
      SimulationData & sim;
  KernelDivergence(SimulationData & s): sim(s){}
  const std::array<int, 3> stencil_start = {-1,-1,-1}, stencil_end = {2, 2, 2};
  const cubism::StencilInfo stencil{-1,-1,-1, 2,2,2, false, {FE_U,FE_V,FE_W,FE_TMPU}};

  template <typename Lab, typename BlockType>
  void operator()(Lab & lab, const cubism::BlockInfo& info, BlockType& o) const
  {
    const Real fac=0.5*info.h*info.h;
    for (int iz=0; iz<FluidBlock::sizeZ; ++iz)
    for (int iy=0; iy<FluidBlock::sizeY; ++iy)
    for (int ix=0; ix<FluidBlock::sizeX; ++ix)
    {
      o(ix,iy,iz).tmpU = (1.0 - o(ix,iy,iz).chi)*fac * ( lab(ix+1,iy,iz).u - lab(ix-1,iy,iz).u +
                                                         lab(ix,iy+1,iz).v - lab(ix,iy-1,iz).v +
                                                         lab(ix,iy,iz+1).w - lab(ix,iy,iz-1).w );
    }

    BlockCase<BlockType> * tempCase = (BlockCase<BlockType> *)(info.auxiliary);
    typename BlockType::ElementType * faceXm = nullptr;
    typename BlockType::ElementType * faceXp = nullptr;
    typename BlockType::ElementType * faceYm = nullptr;
    typename BlockType::ElementType * faceYp = nullptr;
    typename BlockType::ElementType * faceZp = nullptr;
    typename BlockType::ElementType * faceZm = nullptr;
    if (tempCase != nullptr)
    {
      faceXm = tempCase -> storedFace[0] ?  & tempCase -> m_pData[0][0] : nullptr;
      faceXp = tempCase -> storedFace[1] ?  & tempCase -> m_pData[1][0] : nullptr;
      faceYm = tempCase -> storedFace[2] ?  & tempCase -> m_pData[2][0] : nullptr;
      faceYp = tempCase -> storedFace[3] ?  & tempCase -> m_pData[3][0] : nullptr;
      faceZm = tempCase -> storedFace[4] ?  & tempCase -> m_pData[4][0] : nullptr;
      faceZp = tempCase -> storedFace[5] ?  & tempCase -> m_pData[5][0] : nullptr;
    }
    if (faceXm != nullptr)
    {
      int ix = 0;
      for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
      for(int iy=0; iy<FluidBlock::sizeY; ++iy)
      {
        faceXm[iy + FluidBlock::sizeY * iz].clear();
        faceXm[iy + FluidBlock::sizeY * iz].tmpU = (1.0 - o(ix,iy,iz).chi)*fac *(lab(ix-1,iy,iz).u + lab(ix,iy,iz).u);
      }
    }
    if (faceXp != nullptr)
    {
      int ix = FluidBlock::sizeX-1;
      for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
      for(int iy=0; iy<FluidBlock::sizeY; ++iy)
      {
        faceXp[iy + FluidBlock::sizeY * iz].clear();
        faceXp[iy + FluidBlock::sizeY * iz].tmpU = - (1.0 - o(ix,iy,iz).chi)*fac *(lab(ix+1,iy,iz).u + lab(ix,iy,iz).u);
      }
    }
    if (faceYm != nullptr)
    {
      int iy = 0;
      for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
      for(int ix=0; ix<FluidBlock::sizeX; ++ix)
      {
        faceYm[ix + FluidBlock::sizeX * iz].clear();
        faceYm[ix + FluidBlock::sizeX * iz].tmpU = (1.0 - o(ix,iy,iz).chi)*fac *(lab(ix,iy-1,iz).v + lab(ix,iy,iz).v);
      }
    }
    if (faceYp != nullptr)
    {
      int iy = FluidBlock::sizeY-1;
      for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
      for(int ix=0; ix<FluidBlock::sizeX; ++ix)
      {
        faceYp[ix + FluidBlock::sizeX * iz].clear();
        faceYp[ix + FluidBlock::sizeX * iz].tmpU = - (1.0 - o(ix,iy,iz).chi)*fac *(lab(ix,iy+1,iz).v + lab(ix,iy,iz).v);
      }
    }
    if (faceZm != nullptr)
    {
      int iz = 0;
      for(int iy=0; iy<FluidBlock::sizeY; ++iy)
      for(int ix=0; ix<FluidBlock::sizeX; ++ix)
      {
        faceZm[ix + FluidBlock::sizeX * iy].clear();
        faceZm[ix + FluidBlock::sizeX * iy].tmpU = (1.0 - o(ix,iy,iz).chi)*fac *(lab(ix,iy,iz-1).w + lab(ix,iy,iz).w);
      }
    }
    if (faceZp != nullptr)
    {
      int iz = FluidBlock::sizeZ-1;
      for(int iy=0; iy<FluidBlock::sizeY; ++iy)
      for(int ix=0; ix<FluidBlock::sizeX; ++ix)
      {
        faceZp[ix + FluidBlock::sizeX * iy].clear();
        faceZp[ix + FluidBlock::sizeX * iy].tmpU = - (1.0 - o(ix,iy,iz).chi)*fac *(lab(ix,iy,iz+1).w + lab(ix,iy,iz).w);
      }
    }
  }
};

class ComputeDivergence : public Operator
{
  public:
  ComputeDivergence(SimulationData & s) : Operator(s) { }
  void operator()(const double dt)
  {
    sim.startProfiler("Divergence Kernel");

    const KernelDivergence K(sim);
    compute<KernelDivergence>(K,true);

    double div_loc = 0.0;
    const std::vector<cubism::BlockInfo>& myInfo = sim.vInfo();
    #pragma omp parallel for schedule(static) reduction(+: div_loc)
    for(size_t i=0; i<myInfo.size(); i++)
    {
      const cubism::BlockInfo& info = myInfo[i];
      const FluidBlock& b = *(const FluidBlock *)info.ptrBlock;
      for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
      for(int iy=0; iy<FluidBlock::sizeY; ++iy)
      for(int ix=0; ix<FluidBlock::sizeX; ++ix)
        div_loc += std::fabs(b(ix,iy,iz).tmpU);
    }
    double div_tot = 0.0;
    MPI_Reduce(&div_loc, &div_tot, 1, MPI_DOUBLE, MPI_SUM, 0, sim.app_comm);

    size_t loc = myInfo.size();
    size_t tot;
    MPI_Reduce(&loc, &tot, 1, MPI_LONG, MPI_SUM, 0, sim.app_comm);
    if (sim.rank == 0)
    {
      std::cout << "Total div = " << div_tot << std::endl;
      std::ofstream outfile;
      outfile.open("div.txt", std::ios_base::app);
      outfile << sim.time << " " << div_tot << " " << tot<< "\n";
      outfile.close();
    }
    sim.stopProfiler();
  }
  std::string getName() { return "Divergence"; }
};

CubismUP_3D_NAMESPACE_END
#endif // CubismUP_3D_ProcessOperators_h
