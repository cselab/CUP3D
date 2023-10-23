//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//

#pragma once

#include "../SimulationData.h"
#include "Operator.h"

using namespace cubism;

CubismUP_3D_NAMESPACE_BEGIN

//used for mesh refinement
struct GradChiOnTmp
{
  GradChiOnTmp(const SimulationData & s) : sim(s) {}
  const SimulationData & sim;
  const StencilInfo stencil{-2, -2, -2, 3, 3, 3, true, {0}};
  void operator()(ScalarLab & lab, const BlockInfo& info) const
  {
    auto& __restrict__ TMP = (*sim.tmpV)(info.blockID);
    if (info.level == sim.levelMaxVorticity-1 && sim.levelMaxVorticity < sim.levelMax) //allow for refinement up to levelMaxVorticity for vorticity (and up to levelMax-1 for grad chi)
      for(int z=0; z<VectorBlock::sizeZ; ++z)
      for(int y=0; y<VectorBlock::sizeY; ++y)
      for(int x=0; x<VectorBlock::sizeX; ++x)
      {
        if (TMP(x,y,z).magnitude() >= sim.Rtol)
        {
           //set to a value that will not cause refinement or compression
           TMP(x,y,z).u[0] = 0.5*(sim.Rtol+sim.Ctol);
           TMP(x,y,z).u[1] = 0.0;
           TMP(x,y,z).u[2] = 0.0;
        }
      }

    bool done = false;
    const int offset = (info.level == sim.chi->getlevelMax()-1) ? 2 : 1;
    for(int z=-offset; z<VectorBlock::sizeZ+offset; ++z)
    for(int y=-offset; y<VectorBlock::sizeY+offset; ++y)
    for(int x=-offset; x<VectorBlock::sizeX+offset; ++x)
    {
      if (done) break;
      lab(x,y,z).s = std::min(lab(x,y,z).s,(Real)1.0);
      lab(x,y,z).s = std::max(lab(x,y,z).s,(Real)0.0);
      if (lab(x,y,z).s > 0.00001 && lab(x,y,z).s < 0.9)
      {
        TMP(VectorBlock::sizeX/2-1,VectorBlock::sizeY/2-1,VectorBlock::sizeZ/2-1).u[0] = 1e10;
        TMP(VectorBlock::sizeX/2  ,VectorBlock::sizeY/2-1,VectorBlock::sizeZ/2-1).u[0] = 1e10;
        TMP(VectorBlock::sizeX/2-1,VectorBlock::sizeY/2  ,VectorBlock::sizeZ/2-1).u[0] = 1e10;
        TMP(VectorBlock::sizeX/2  ,VectorBlock::sizeY/2-1,VectorBlock::sizeZ/2-1).u[0] = 1e10;
        TMP(VectorBlock::sizeX/2-1,VectorBlock::sizeY/2-1,VectorBlock::sizeZ/2  ).u[0] = 1e10;
        TMP(VectorBlock::sizeX/2  ,VectorBlock::sizeY/2  ,VectorBlock::sizeZ/2  ).u[0] = 1e10;
        TMP(VectorBlock::sizeX/2-1,VectorBlock::sizeY/2-1,VectorBlock::sizeZ/2  ).u[0] = 1e10;
        TMP(VectorBlock::sizeX/2  ,VectorBlock::sizeY/2-1,VectorBlock::sizeZ/2  ).u[0] = 1e10;
        done = true;
        break;
      }
      else if (lab(x,y,z).s > 0.9 && z >= 0 && z <VectorBlock::sizeZ && y >= 0 && y <VectorBlock::sizeY && x >= 0 && x <VectorBlock::sizeX) //compress the grid if inside an obstacle
      {
        TMP(x,y,z).u[0] = 0.0;
        TMP(x,y,z).u[1] = 0.0;
        TMP(x,y,z).u[2] = 0.0;
      }
    }
  }
};

inline Real findMaxU(SimulationData& sim)
{
  const std::vector<BlockInfo>& myInfo = sim.velInfo();
  const Real uinf[3] = {sim.uinf[0], sim.uinf[1], sim.uinf[2]};
  Real maxU = 0;
  #pragma omp parallel for schedule(static) reduction(max : maxU)
  for(size_t i=0; i<myInfo.size(); i++)
  {
    const VectorBlock& b = *(const VectorBlock *)myInfo[i].ptrBlock;
    for(int z=0; z<VectorBlock::sizeZ; ++z)
    for(int y=0; y<VectorBlock::sizeY; ++y)
    for(int x=0; x<VectorBlock::sizeX; ++x)
    {
      const Real advu = std::fabs(b(x,y,z).u[0] + uinf[0]);
      const Real advv = std::fabs(b(x,y,z).u[1] + uinf[1]);
      const Real advw = std::fabs(b(x,y,z).u[2] + uinf[2]);
      const Real maxUl = std::max({advu, advv, advw});
      maxU = std::max(maxU, maxUl);
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, & maxU, 1, MPI_Real, MPI_MAX, sim.comm);
  assert(maxU >= 0);
  return maxU;
}

struct KernelVorticity
{
  SimulationData & sim;
  KernelVorticity(SimulationData & s): sim(s){};
  const StencilInfo stencil{-1,-1,-1, 2,2,2, false, {0,1,2}};
  const std::vector<BlockInfo> & vInfo = sim.tmpVInfo();
  const int Nx = VectorBlock::sizeX;
  const int Ny = VectorBlock::sizeY;
  const int Nz = VectorBlock::sizeZ;

  void operator()(const VectorLab & lab, const BlockInfo& info) const
  {
    const cubism::BlockInfo& info2 = vInfo[info.blockID];
    VectorBlock& o = *(VectorBlock*)info2.ptrBlock;
    const Real inv2h = .5 * info.h * info.h;
    for (int z=0; z<Nz; ++z)
    for (int y=0; y<Ny; ++y)
    for (int x=0; x<Nx; ++x)
    {
      const VectorElement &LW=lab(x-1,y,z), &LE=lab(x+1,y,z);
      const VectorElement &LS=lab(x,y-1,z), &LN=lab(x,y+1,z);
      const VectorElement &LF=lab(x,y,z-1), &LB=lab(x,y,z+1);
      o(x,y,z).u[0] = inv2h * ( (LN.u[2]-LS.u[2]) - (LB.u[1]-LF.u[1]) );
      o(x,y,z).u[1] = inv2h * ( (LB.u[0]-LF.u[0]) - (LE.u[2]-LW.u[2]) );
      o(x,y,z).u[2] = inv2h * ( (LE.u[1]-LW.u[1]) - (LN.u[0]-LS.u[0]) );
    }
    BlockCase<VectorBlock> * tempCase = (BlockCase<VectorBlock> *)(info.auxiliary);

    if (tempCase == nullptr) return;

    VectorElement * const faceXm = tempCase -> storedFace[0] ?  & tempCase -> m_pData[0][0] : nullptr;
    VectorElement * const faceXp = tempCase -> storedFace[1] ?  & tempCase -> m_pData[1][0] : nullptr;
    VectorElement * const faceYm = tempCase -> storedFace[2] ?  & tempCase -> m_pData[2][0] : nullptr;
    VectorElement * const faceYp = tempCase -> storedFace[3] ?  & tempCase -> m_pData[3][0] : nullptr;
    VectorElement * const faceZm = tempCase -> storedFace[4] ?  & tempCase -> m_pData[4][0] : nullptr;
    VectorElement * const faceZp = tempCase -> storedFace[5] ?  & tempCase -> m_pData[5][0] : nullptr;
    if (faceXm != nullptr)
    {
       const int x = 0;
       for(int z=0; z<Nz; ++z)
       for(int y=0; y<Ny; ++y)
       {
          const VectorElement &LW=lab(x-1,y,z);
          const VectorElement &LC=lab(x,y,z);
          faceXm[y + Ny * z].u[1] = -inv2h*( LW.u[2] + LC.u[2] );
          faceXm[y + Ny * z].u[2] = +inv2h*( LW.u[1] + LC.u[1] );
       }
    }
    if (faceXp != nullptr)
    {
       const int x = Nx-1;
       for(int z=0; z<Nz; ++z)
       for(int y=0; y<Ny; ++y)
       {
          const VectorElement &LE=lab(x+1,y,z);
          const VectorElement &LC=lab(x,y,z);
          faceXp[y + Ny * z].u[1] = +inv2h*( LE.u[2] + LC.u[2] );
          faceXp[y + Ny * z].u[2] = -inv2h*( LE.u[1] + LC.u[1] );
       }
    }
    if (faceYm != nullptr)
    {
       const int y = 0;
       for(int z=0; z<Nz; ++z)
       for(int x=0; x<Nx; ++x)
       {
          const VectorElement &LS=lab(x,y-1,z);
          const VectorElement &LC=lab(x,y,z);
          faceYm[x + Nx * z].u[0] = +inv2h*(LS.u[2]+LC.u[2]);
          faceYm[x + Nx * z].u[2] = -inv2h*(LS.u[0]+LC.u[0]);
       }
     }
     if (faceYp != nullptr)
     {
       const int y = Ny-1;
       for(int z=0; z<Nz; ++z)
       for(int x=0; x<Nx; ++x)
       {
          const VectorElement &LN=lab(x,y+1,z);
          const VectorElement &LC=lab(x,y,z);
          faceYp[x + Nx * z].u[0] = -inv2h*(LN.u[2]+LC.u[2]);
          faceYp[x + Nx * z].u[2] = +inv2h*(LN.u[0]+LC.u[0]);
       }
     }
     if (faceZm != nullptr)
     {
       const int z = 0;
       for(int y=0; y<Ny; ++y)
       for(int x=0; x<Nx; ++x)
       {
          const VectorElement &LF=lab(x,y,z-1);
          const VectorElement &LC=lab(x,y,z);
          faceZm[x + Nx * y].u[0] = -inv2h*(LF.u[1]+LC.u[1]);
          faceZm[x + Nx * y].u[1] = +inv2h*(LF.u[0]+LC.u[0]);
       }
     }
     if (faceZp != nullptr)
     {
       const int z = Nz-1;
       for(int y=0; y<Ny; ++y)
       for(int x=0; x<Nx; ++x)
       {
          const VectorElement &LB=lab(x,y,z+1);
          const VectorElement &LC=lab(x,y,z);
          faceZp[x + Nx * y].u[0] = +inv2h*(LB.u[1]+LC.u[1]);
          faceZp[x + Nx * y].u[1] = -inv2h*(LB.u[0]+LC.u[0]);
       }
     }
  }
};

class ComputeVorticity : public Operator
{
  public:
  ComputeVorticity(SimulationData & s) : Operator(s) { }
  void operator()(const Real dt)
  {
    const KernelVorticity K(sim);
    compute<VectorLab>(K,sim.vel,sim.tmpV);
    const std::vector<BlockInfo>& myInfo = sim.tmpVInfo();
    #ifdef PRESERVE_SYMMETRY
    Real omega_max_min [6] = {0.};
    #pragma omp parallel for reduction (max:omega_max_min[:6])
    #else
    #pragma omp parallel for
    #endif
    for(size_t i=0; i<myInfo.size(); i++)
    {
      const BlockInfo& info = myInfo[i];
      VectorBlock& b = *( VectorBlock *)info.ptrBlock;
      const Real fac = 1.0/(info.h*info.h*info.h);
      for(int z=0; z<VectorBlock::sizeZ; ++z)
      for(int y=0; y<VectorBlock::sizeY; ++y)
      for(int x=0; x<VectorBlock::sizeX; ++x)
      {
        b(x,y,z).u[0] *=fac;
        b(x,y,z).u[1] *=fac;
        b(x,y,z).u[2] *=fac;
        #ifdef PRESERVE_SYMMETRY
        omega_max_min[0] = std::max(omega_max_min[0], b(x,y,z).u[0]);
        omega_max_min[1] = std::max(omega_max_min[1], b(x,y,z).u[1]);
        omega_max_min[2] = std::max(omega_max_min[2], b(x,y,z).u[2]);
        omega_max_min[3] = std::max(omega_max_min[3],-b(x,y,z).u[0]);
        omega_max_min[4] = std::max(omega_max_min[4],-b(x,y,z).u[1]);
        omega_max_min[5] = std::max(omega_max_min[5],-b(x,y,z).u[2]);
        #endif
      }
    }
    #ifdef PRESERVE_SYMMETRY
    MPI_Reduce(sim.rank == 0 ? MPI_IN_PLACE: omega_max_min, omega_max_min, 6, MPI_Real, MPI_MAX, 0, sim.comm);
    if (sim.rank == 0 && sim.verbose)
    {
       std::cout << "Vorticity (x): max=" << omega_max_min[0] << " min=" << -omega_max_min[3] << " difference: " << omega_max_min[0]-omega_max_min[3] << std::endl;
       std::cout << "Vorticity (y): max=" << omega_max_min[1] << " min=" << -omega_max_min[4] << " difference: " << omega_max_min[1]-omega_max_min[4] << std::endl;
       std::cout << "Vorticity (z): max=" << omega_max_min[2] << " min=" << -omega_max_min[5] << " difference: " << omega_max_min[2]-omega_max_min[5] << std::endl;
    }
    #endif
  }
  std::string getName() { return "Vorticity"; }
};

class KernelQcriterion
{
  public:
  SimulationData & sim;
  KernelQcriterion(SimulationData & s): sim(s){};
  const std::array<int, 3> stencil_start = {-1,-1,-1}, stencil_end = {2, 2, 2};
  const cubism::StencilInfo stencil{-1,-1,-1, 2,2,2, false, {0,1,2}};
  const std::vector<cubism::BlockInfo> & vInfo = sim.presInfo();

  void operator()(VectorLab & lab, const cubism::BlockInfo& info) const
  {
    ScalarBlock& o = *( ScalarBlock *)vInfo[info.blockID].ptrBlock;
    const Real inv2h = .5 / info.h;
    for (int iz=0; iz<ScalarBlock::sizeZ; ++iz)
    for (int iy=0; iy<ScalarBlock::sizeY; ++iy)
    for (int ix=0; ix<ScalarBlock::sizeX; ++ix) {
      const VectorElement &LW=lab(ix-1,iy,iz), &LE=lab(ix+1,iy,iz);
      const VectorElement &LS=lab(ix,iy-1,iz), &LN=lab(ix,iy+1,iz);
      const VectorElement &LF=lab(ix,iy,iz-1), &LB=lab(ix,iy,iz+1);
      const Real WX  = inv2h * ( (LN.u[2]-LS.u[2]) - (LB.u[1]-LF.u[1]) );
      const Real WY  = inv2h * ( (LB.u[0]-LF.u[0]) - (LE.u[2]-LW.u[2]) );
      const Real WZ  = inv2h * ( (LE.u[1]-LW.u[1]) - (LN.u[0]-LS.u[0]) );
      const Real D11 = inv2h * (LE.u[0]-LW.u[0]); // shear stresses
      const Real D22 = inv2h * (LN.u[1]-LS.u[1]); // shear stresses
      const Real D33 = inv2h * (LB.u[2]-LF.u[2]); // shear stresses
      const Real D12 = inv2h * (LN.u[0]-LS.u[0] + LE.u[1]-LW.u[1]); // shear stresses
      const Real D13 = inv2h * (LE.u[2]-LW.u[2] + LB.u[0]-LF.u[0]); // shear stresses
      const Real D23 = inv2h * (LB.u[1]-LF.u[1] + LN.u[2]-LS.u[2]); // shear stresses
      // trace( S S^t ) where S is the sym part of the vel gradient:
      const Real SS = D11*D11 +D22*D22 +D33*D33 +(D12*D12 +D13*D13 +D23*D23)/2;
      o(ix,iy,iz).s = ( (WX*WX + WY*WY + WZ*WZ)/2 - SS ) / 2;
    }
  }
};

class ComputeQcriterion : public Operator
{
  public:
  ComputeQcriterion(SimulationData & s) : Operator(s) { }
  void operator()(const Real dt)
  {
    const KernelQcriterion K(sim);
    cubism::compute<VectorLab>(K,sim.vel);
  }
  std::string getName() { return "Qcriterion"; }
};

class KernelDivergence
{
  public:
  SimulationData & sim;
  KernelDivergence(SimulationData & s): sim(s){}
  const std::array<int, 3> stencil_start = {-1,-1,-1}, stencil_end = {2, 2, 2};
  const cubism::StencilInfo stencil{-1,-1,-1, 2,2,2, false, {0,1,2}};
  const std::vector<cubism::BlockInfo> & vInfo = sim.tmpVInfo();
  const std::vector<cubism::BlockInfo> & chiInfo = sim.chiInfo();

  void operator()(VectorLab & lab, const cubism::BlockInfo& info) const
  {
    VectorBlock& o = *( VectorBlock *)vInfo[info.blockID].ptrBlock;
    ScalarBlock& c = *( ScalarBlock *)chiInfo[info.blockID].ptrBlock;
    const Real fac=0.5*info.h*info.h;
    for (int iz=0; iz<VectorBlock::sizeZ; ++iz)
    for (int iy=0; iy<VectorBlock::sizeY; ++iy)
    for (int ix=0; ix<VectorBlock::sizeX; ++ix)
    {
      o(ix,iy,iz).u[0] = (1.0 - c(ix,iy,iz).s)*fac * ( lab(ix+1,iy,iz).u[0] - lab(ix-1,iy,iz).u[0] +
                                                         lab(ix,iy+1,iz).u[1] - lab(ix,iy-1,iz).u[1] +
                                                         lab(ix,iy,iz+1).u[2] - lab(ix,iy,iz-1).u[2] );
    }

    BlockCase<VectorBlock> * tempCase = (BlockCase<VectorBlock> *)(info.auxiliary);
    typename VectorBlock::ElementType * faceXm = nullptr;
    typename VectorBlock::ElementType * faceXp = nullptr;
    typename VectorBlock::ElementType * faceYm = nullptr;
    typename VectorBlock::ElementType * faceYp = nullptr;
    typename VectorBlock::ElementType * faceZp = nullptr;
    typename VectorBlock::ElementType * faceZm = nullptr;
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
      for(int iz=0; iz<ScalarBlock::sizeZ; ++iz)
      for(int iy=0; iy<ScalarBlock::sizeY; ++iy)
      {
        faceXm[iy + ScalarBlock::sizeY * iz].clear();
        faceXm[iy + ScalarBlock::sizeY * iz].u[0] = (1.0 - c(ix,iy,iz).s)*fac *(lab(ix-1,iy,iz).u[0] + lab(ix,iy,iz).u[0]);
      }
    }
    if (faceXp != nullptr)
    {
      int ix = ScalarBlock::sizeX-1;
      for(int iz=0; iz<ScalarBlock::sizeZ; ++iz)
      for(int iy=0; iy<ScalarBlock::sizeY; ++iy)
      {
        faceXp[iy + ScalarBlock::sizeY * iz].clear();
        faceXp[iy + ScalarBlock::sizeY * iz].u[0] = - (1.0 - c(ix,iy,iz).s)*fac *(lab(ix+1,iy,iz).u[0] + lab(ix,iy,iz).u[0]);
      }
    }
    if (faceYm != nullptr)
    {
      int iy = 0;
      for(int iz=0; iz<ScalarBlock::sizeZ; ++iz)
      for(int ix=0; ix<ScalarBlock::sizeX; ++ix)
      {
        faceYm[ix + ScalarBlock::sizeX * iz].clear();
        faceYm[ix + ScalarBlock::sizeX * iz].u[0] = (1.0 - c(ix,iy,iz).s)*fac *(lab(ix,iy-1,iz).u[1] + lab(ix,iy,iz).u[1]);
      }
    }
    if (faceYp != nullptr)
    {
      int iy = ScalarBlock::sizeY-1;
      for(int iz=0; iz<ScalarBlock::sizeZ; ++iz)
      for(int ix=0; ix<ScalarBlock::sizeX; ++ix)
      {
        faceYp[ix + ScalarBlock::sizeX * iz].clear();
        faceYp[ix + ScalarBlock::sizeX * iz].u[0] = - (1.0 - c(ix,iy,iz).s)*fac *(lab(ix,iy+1,iz).u[1] + lab(ix,iy,iz).u[1]);
      }
    }
    if (faceZm != nullptr)
    {
      int iz = 0;
      for(int iy=0; iy<ScalarBlock::sizeY; ++iy)
      for(int ix=0; ix<ScalarBlock::sizeX; ++ix)
      {
        faceZm[ix + ScalarBlock::sizeX * iy].clear();
        faceZm[ix + ScalarBlock::sizeX * iy].u[0] = (1.0 - c(ix,iy,iz).s)*fac *(lab(ix,iy,iz-1).u[2] + lab(ix,iy,iz).u[2]);
      }
    }
    if (faceZp != nullptr)
    {
      int iz = ScalarBlock::sizeZ-1;
      for(int iy=0; iy<ScalarBlock::sizeY; ++iy)
      for(int ix=0; ix<ScalarBlock::sizeX; ++ix)
      {
        faceZp[ix + ScalarBlock::sizeX * iy].clear();
        faceZp[ix + ScalarBlock::sizeX * iy].u[0] = - (1.0 - c(ix,iy,iz).s)*fac *(lab(ix,iy,iz+1).u[2] + lab(ix,iy,iz).u[2]);
      }
    }
  }
};

class ComputeDivergence : public Operator
{
  public:
  ComputeDivergence(SimulationData & s) : Operator(s) { }
  void operator()(const Real dt)
  {
    const KernelDivergence K(sim);
    cubism::compute<VectorLab>(K,sim.vel,sim.chi);

    Real div_loc = 0.0;
    const std::vector<cubism::BlockInfo>& myInfo = sim.tmpVInfo();
    #pragma omp parallel for schedule(static) reduction(+: div_loc)
    for(size_t i=0; i<myInfo.size(); i++)
    {
      const cubism::BlockInfo& info = myInfo[i];
      const VectorBlock& b = *(const VectorBlock *)info.ptrBlock;
      for(int iz=0; iz<VectorBlock::sizeZ; ++iz)
      for(int iy=0; iy<VectorBlock::sizeY; ++iy)
      for(int ix=0; ix<VectorBlock::sizeX; ++ix)
        div_loc += std::fabs(b(ix,iy,iz).u[0]);
    }
    Real div_tot = 0.0;
    MPI_Reduce(&div_loc, &div_tot, 1, MPI_Real, MPI_SUM, 0, sim.comm);

    size_t loc = myInfo.size();
    size_t tot;
    MPI_Reduce(&loc, &tot, 1, MPI_LONG, MPI_SUM, 0, sim.comm);
    if (sim.rank == 0)
    {
      std::cout << "Total div = " << div_tot << std::endl;
      std::ofstream outfile;
      outfile.open("div.txt", std::ios_base::app);
      outfile << sim.time << " " << div_tot << " " << tot<< "\n";
      outfile.close();
    }
  }
  std::string getName() { return "Divergence"; }
};

CubismUP_3D_NAMESPACE_END
