//
//  CubismUP_3D
//  Copyright (c) 2020 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Michalis Chatzimanolakis (michaich@ethz.ch).
//

#include "PressureProjection.h"

CubismUP_3D_NAMESPACE_BEGIN

struct KernelGradP
{
  const StencilInfo stencil{-1,-1,-1,2,2,2,false,{0}};
  SimulationData & sim;
  const std::vector<BlockInfo>& tmpVInfo = sim.tmpVInfo();
  const Real dt = sim.dt;
  const int Nx = VectorBlock::sizeX;
  const int Ny = VectorBlock::sizeY;
  const int Nz = VectorBlock::sizeZ;

  KernelGradP(SimulationData & s): sim(s) {}

  ~KernelGradP() {}

  void operator()(const ScalarLab & lab, const BlockInfo& info) const
  {
    VectorBlock& o = *(VectorBlock*)tmpVInfo[info.blockID].ptrBlock;
    const Real fac = -0.5*dt*info.h*info.h;
    for(int z=0; z<Nz; ++z)
    for(int y=0; y<Ny; ++y)
    for(int x=0; x<Nx; ++x)
    {
      o(x,y,z).u[0] = fac*(lab(x+1,y,z).s-lab(x-1,y,z).s);
      o(x,y,z).u[1] = fac*(lab(x,y+1,z).s-lab(x,y-1,z).s);
      o(x,y,z).u[2] = fac*(lab(x,y,z+1).s-lab(x,y,z-1).s);
    }
    BlockCase<VectorBlock> * tempCase = (BlockCase<VectorBlock> *)(tmpVInfo[info.blockID].auxiliary);

    if (tempCase == nullptr) return; //no flux corrections needed for this block

    VectorElement * faceXm = tempCase -> storedFace[0] ?  & tempCase -> m_pData[0][0] : nullptr;
    VectorElement * faceXp = tempCase -> storedFace[1] ?  & tempCase -> m_pData[1][0] : nullptr;
    VectorElement * faceYm = tempCase -> storedFace[2] ?  & tempCase -> m_pData[2][0] : nullptr;
    VectorElement * faceYp = tempCase -> storedFace[3] ?  & tempCase -> m_pData[3][0] : nullptr;
    VectorElement * faceZm = tempCase -> storedFace[4] ?  & tempCase -> m_pData[4][0] : nullptr;
    VectorElement * faceZp = tempCase -> storedFace[5] ?  & tempCase -> m_pData[5][0] : nullptr;
    if (faceXm != nullptr)
    {
      const int x = 0;
      for(int z=0; z<Nz; ++z)
      for(int y=0; y<Ny; ++y)
        faceXm[y + Ny * z].u[0] =   fac *(lab(x-1,y,z).s + lab(x,y,z).s);
    }
    if (faceXp != nullptr)
    {
      const int x = Nx-1;
      for(int z=0; z<Nz; ++z)
      for(int y=0; y<Ny; ++y)
        faceXp[y + Ny * z].u[0] = - fac *(lab(x+1,y,z).s + lab(x,y,z).s);
    }
    if (faceYm != nullptr)
    {
      const int y = 0;
      for(int z=0; z<Nz; ++z)
      for(int x=0; x<Nx; ++x)
        faceYm[x + Nx * z].u[1] =   fac *(lab(x,y-1,z).s + lab(x,y,z).s);
    }
    if (faceYp != nullptr)
    {
      const int y = Ny-1;
      for(int z=0; z<Nz; ++z)
      for(int x=0; x<Nx; ++x)
        faceYp[x + Nx * z].u[1] = - fac *(lab(x,y+1,z).s + lab(x,y,z).s);
    }
    if (faceZm != nullptr)
    {
      const int z = 0;
      for(int y=0; y<Ny; ++y)
      for(int x=0; x<Nx; ++x)
        faceZm[x + Nx * y].u[2] =   fac *(lab(x,y,z-1).s + lab(x,y,z).s);
    }
    if (faceZp != nullptr)
    {
      const int z = Nz-1;
      for(int y=0; y<Ny; ++y)
      for(int x=0; x<Nx; ++x)
        faceZp[x + Nx * y].u[2] = - fac *(lab(x,y,z+1).s + lab(x,y,z).s);
    }
  }
};

PressureProjection::PressureProjection(SimulationData & s) : Operator(s)
{
  pressureSolver = makePoissonSolver(s);
  sim.pressureSolver = pressureSolver;
}

void PressureProjection::operator()(const Real dt)
{
  //Solve the Poisson equation and use pressure to perform
  //pressure projection of the velocity field.

  //The rhs of the linear system is contained in sim.lhs
  //The initial guess is contained in sim.pres
  //Here we solve for phi := p^{n+1}-p^{n} if step < step_2nd_start
  pressureSolver->solve();

  const int Nx = VectorBlock::sizeX;
  const int Ny = VectorBlock::sizeY;
  const int Nz = VectorBlock::sizeZ;
  const std::vector<BlockInfo>& velInfo  = sim.velInfo();
  const std::vector<BlockInfo>& tmpVInfo = sim.tmpVInfo();
  const std::vector<BlockInfo>& presInfo = sim.presInfo();
  const std::vector<BlockInfo>& pOldInfo = sim.pOldInfo();

  if (sim.step > sim.step_2nd_start) //recover p^{n+1} = phi + p^{n}
  {
    #pragma omp parallel for
    for(size_t i=0; i<presInfo.size(); i++)
    {
      ScalarBlock& p          = *(ScalarBlock*)presInfo[i].ptrBlock;
      const ScalarBlock& pOld = *(ScalarBlock*)pOldInfo[i].ptrBlock;
      for (int z=0; z<Nz; ++z)
      for (int y=0; y<Ny; ++y)
      for (int x=0; x<Nx; ++x)
        p(x,y,z) += pOld(x,y,z);
    }
  }

  //Compute grad(P) and put it to the vector tmpV
  compute<ScalarLab>(KernelGradP(sim),sim.pres,sim.tmpV);

  //Perform the projection and set u^{n+1} = u - grad(P)*dt
  #pragma omp parallel for
  for(size_t i=0; i<velInfo.size(); i++)
  {
    const Real fac = 1.0/(velInfo[i].h*velInfo[i].h*velInfo[i].h);
    const VectorBlock& gradP = *(VectorBlock*)tmpVInfo[i].ptrBlock;
    VectorBlock& v = *(VectorBlock*)velInfo[i].ptrBlock;
    for (int z=0; z<Nz; ++z)
    for (int y=0; y<Ny; ++y)
    for (int x=0; x<Nx; ++x)
    {
      v(x,y,z).u[0] += fac*gradP(x,y,z).u[0];
      v(x,y,z).u[1] += fac*gradP(x,y,z).u[1];
      v(x,y,z).u[2] += fac*gradP(x,y,z).u[2];
    }
  }
}

CubismUP_3D_NAMESPACE_END
