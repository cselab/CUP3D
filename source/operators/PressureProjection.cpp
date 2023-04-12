//
//  CubismUP_3D
//  Copyright (c) 2020 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Michalis Chatzimanolakis (michaich@ethz.ch).
//

#include "PressureProjection.h"

CubismUP_3D_NAMESPACE_BEGIN

using CHIMAT =  Real[CUP_BLOCK_SIZEZ][CUP_BLOCK_SIZEY][CUP_BLOCK_SIZEX];
using UDEFMAT = Real[CUP_BLOCK_SIZEZ][CUP_BLOCK_SIZEY][CUP_BLOCK_SIZEX][3];

struct KernelDivPressure
{
  const SimulationData & sim;
  const StencilInfo stencil = StencilInfo(-1,-1,-1,2,2,2,false,{0});
  const std::vector<BlockInfo>& tmpVInfo = sim.tmpVInfo();
  const int Nx = VectorBlock::sizeX;
  const int Ny = VectorBlock::sizeY;
  const int Nz = VectorBlock::sizeZ;

  KernelDivPressure(const SimulationData& s) :sim(s) {}

  void operator()(const ScalarLab & lab, const BlockInfo& info) const 
  {
    VectorBlock & __restrict__ b = (*sim.tmpV)(info.blockID);
    const Real fac = info.h;
    for(int z=0; z<Nz; ++z)
    for(int y=0; y<Ny; ++y)
    for(int x=0; x<Nx; ++x)
        #ifdef PRESERVE_SYMMETRY
              b(x,y,z).u[0] = fac*( ConsistentSum(lab(x+1,y,z).s+lab(x-1,y,z).s,lab(x,y+1,z).s+lab(x,y-1,z).s,lab(x,y,z+1).s+lab(x,y,z-1).s) - 6.0*lab(x,y,z).s);
        #else
              b(x,y,z).u[0] = fac*(lab(x+1,y,z).s+lab(x-1,y,z).s +lab(x,y+1,z).s+lab(x,y-1,z).s +lab(x,y,z+1).s+lab(x,y,z-1).s - 6.0*lab(x,y,z).s);
        #endif

    BlockCase<VectorBlock> * tempCase = (BlockCase<VectorBlock> *)(tmpVInfo[info.blockID].auxiliary);

    if (tempCase == nullptr) return; //no flux corrections needed for this block

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
        faceXm[y + Ny * z].u[0] = fac *(lab(x,y,z).s - lab(x-1,y,z).s);
    }
    if (faceXp != nullptr)
    {
      const int x = Nx-1;
      for(int z=0; z<Nz; ++z)
      for(int y=0; y<Ny; ++y)
        faceXp[y + Ny * z].u[0] = - fac *(lab(x+1,y,z).s - lab(x,y,z).s);
    }
    if (faceYm != nullptr)
    {
      const int y = 0;
      for(int z=0; z<Nz; ++z)
      for(int x=0; x<Nx; ++x)
        faceYm[x + Nx * z].u[0] = fac *(lab(x,y,z).s - lab(x,y-1,z).s);
    }
    if (faceYp != nullptr)
    {
      const int y = Ny-1;
      for(int z=0; z<Nz; ++z)
      for(int x=0; x<Nx; ++x)
        faceYp[x + Nx * z].u[0] = - fac *(lab(x,y+1,z).s - lab(x,y,z).s);
    }
    if (faceZm != nullptr)
    {
      const int z = 0;
      for(int y=0; y<Ny; ++y)
      for(int x=0; x<Nx; ++x)
        faceZm[x + Nx * y].u[0] = fac *(lab(x,y,z).s - lab(x,y,z-1).s);
    }
    if (faceZp != nullptr)
    {
      const int z = Nz-1;
      for(int y=0; y<Ny; ++y)
      for(int x=0; x<Nx; ++x)
        faceZp[x + Nx * y].u[0] = - fac *(lab(x,y,z+1).s - lab(x,y,z).s);
    }
  }
};

struct KernelPressureRHS
{
  SimulationData & sim;
  const Real dt;
  ObstacleVector * const obstacle_vector = sim.obstacle_vector;
  const int nShapes = obstacle_vector->nObstacles();
  StencilInfo stencil  = StencilInfo(-1,-1,-1, 2,2,2, false, {0,1,2});
  StencilInfo stencil2 = StencilInfo(-1,-1,-1, 2,2,2, false, {0,1,2});
  const std::vector<BlockInfo>& lhsInfo = sim.lhsInfo();
  const std::vector<BlockInfo>& chiInfo = sim.chiInfo();
  const int Nx = VectorBlock::sizeX;
  const int Ny = VectorBlock::sizeY;
  const int Nz = VectorBlock::sizeZ;

  KernelPressureRHS(SimulationData& s, const Real a_dt) :sim(s), dt(a_dt) {}

  void operator()(const VectorLab & lab, const VectorLab & uDefLab, const BlockInfo& info, const BlockInfo& info2) const
  {
    const Real h = info.h, fac = 0.5*h*h/dt;
    const ScalarBlock & __restrict__ c  = (*sim.chi)(info2.blockID);
    ScalarBlock & __restrict__ p  = (*sim.lhs)(info2.blockID);

    for(int z=0; z<Nz; ++z)
    for(int y=0; y<Ny; ++y)
    for(int x=0; x<Nx; ++x)
    {
      {
        const VectorElement &LW = lab(x-1,y,  z  ), &LE = lab(x+1,y,  z  );
        const VectorElement &LS = lab(x,  y-1,z  ), &LN = lab(x,  y+1,z  );
        const VectorElement &LF = lab(x,  y,  z-1), &LB = lab(x,  y,  z+1);
        #ifdef PRESERVE_SYMMETRY
        p(x,y,z).s = fac*ConsistentSum(LE.u[0]-LW.u[0], LN.u[1]-LS.u[1], LB.u[2]-LF.u[2]);
        #else
        p(x,y,z).s = fac*(LE.u[0]-LW.u[0] + LN.u[1]-LS.u[1] + LB.u[2]-LF.u[2]);
        #endif
      }
      {
        const VectorElement &LW = uDefLab(x-1,y,  z  ), &LE = uDefLab(x+1,y,  z  );
        const VectorElement &LS = uDefLab(x,  y-1,z  ), &LN = uDefLab(x,  y+1,z  );
        const VectorElement &LF = uDefLab(x,  y,  z-1), &LB = uDefLab(x,  y,  z+1);
        #ifdef PRESERVE_SYMMETRY
        const Real divUs = ConsistentSum(LE.u[0]-LW.u[0], LN.u[1]-LS.u[1], LB.u[2]-LF.u[2]);
        #else
        const Real divUs = LE.u[0]-LW.u[0] + LN.u[1]-LS.u[1] + LB.u[2]-LF.u[2];
        #endif
        p(x,y,z).s += - c(x,y,z).s * fac * divUs;
      }
    }

    BlockCase<ScalarBlock> * tempCase = (BlockCase<ScalarBlock> *)(lhsInfo[info2.blockID].auxiliary);

    if (tempCase == nullptr) return; //no flux corrections needed for this block

    ScalarElement * const faceXm = tempCase -> storedFace[0] ?  & tempCase -> m_pData[0][0] : nullptr;
    ScalarElement * const faceXp = tempCase -> storedFace[1] ?  & tempCase -> m_pData[1][0] : nullptr;
    ScalarElement * const faceYm = tempCase -> storedFace[2] ?  & tempCase -> m_pData[2][0] : nullptr;
    ScalarElement * const faceYp = tempCase -> storedFace[3] ?  & tempCase -> m_pData[3][0] : nullptr;
    ScalarElement * const faceZm = tempCase -> storedFace[4] ?  & tempCase -> m_pData[4][0] : nullptr;
    ScalarElement * const faceZp = tempCase -> storedFace[5] ?  & tempCase -> m_pData[5][0] : nullptr;

    if (faceXm != nullptr)
    {
      const int x = 0;
      for(int z=0; z<Nz; ++z)
      for(int y=0; y<Ny; ++y)
        faceXm[y + Ny * z].s =  fac*(lab(x-1,y,z).u[0] + lab(x,y,z).u[0]) - c(x,y,z).s*fac*(uDefLab(x-1,y,z).u[0] + uDefLab(x,y,z).u[0]);
    }
    if (faceXp != nullptr)
    {
      const int x = Nx-1;
      for(int z=0; z<Nz; ++z)
      for(int y=0; y<Ny; ++y)
        faceXp[y + Ny * z].s = -fac*(lab(x+1,y,z).u[0] + lab(x,y,z).u[0]) + c(x,y,z).s*fac*(uDefLab(x+1,y,z).u[0] + uDefLab(x,y,z).u[0]);
    }
    if (faceYm != nullptr)
    {
      const int y = 0;
      for(int z=0; z<Nz; ++z)
      for(int x=0; x<Nx; ++x)
        faceYm[x + Nx * z].s =  fac*(lab(x,y-1,z).u[1] + lab(x,y,z).u[1]) - c(x,y,z).s*fac*(uDefLab(x,y-1,z).u[1] + uDefLab(x,y,z).u[1]);
    }
    if (faceYp != nullptr)
    {
      const int y = Ny-1;
      for(int z=0; z<Nz; ++z)
      for(int x=0; x<Nx; ++x)
        faceYp[x + Nx * z].s = -fac*(lab(x,y+1,z).u[1] + lab(x,y,z).u[1]) + c(x,y,z).s*fac*(uDefLab(x,y+1,z).u[1] + uDefLab(x,y,z).u[1]);
    }
    if (faceZm != nullptr)
    {
      const int z = 0;
      for(int y=0; y<Ny; ++y)
      for(int x=0; x<Nx; ++x)
        faceZm[x + Nx * y].s =  fac*(lab(x,y,z-1).u[2] + lab(x,y,z).u[2]) - c(x,y,z).s*fac*(uDefLab(x,y,z-1).u[2] + uDefLab(x,y,z).u[2]);
    }
    if (faceZp != nullptr)
    {
      const int z = Nz-1;
      for(int y=0; y<Ny; ++y)
      for(int x=0; x<Nx; ++x)
        faceZp[x + Nx * y].s = -fac*(lab(x,y,z+1).u[2] + lab(x,y,z).u[2]) + c(x,y,z).s*fac*(uDefLab(x,y,z+1).u[2] + uDefLab(x,y,z).u[2]);
    }
  }
};

/// Add obstacle's udef to tmpV
static void kernelUpdateTmpV(SimulationData& sim)
{
  const int Nx = VectorBlock::sizeX;
  const int Ny = VectorBlock::sizeY;
  const int Nz = VectorBlock::sizeZ;
  const std::vector<BlockInfo>& chiInfo = sim.chiInfo();
  #pragma omp parallel
  {
    for (const auto &obstacle : sim.obstacle_vector->getObstacleVector())
    {
      const auto& obstblocks = obstacle->getObstacleBlocks();
      #pragma omp for schedule(dynamic, 1)
      for (size_t i = 0; i < chiInfo.size(); ++i)
      {
        const BlockInfo& info = chiInfo[i];
        const auto pos = obstblocks[info.blockID];
        if(pos == nullptr) continue;

        const ScalarBlock& c = (*sim.chi)(i);
        VectorBlock& b = (*sim.tmpV)(i);
        const UDEFMAT & __restrict__ UDEF = pos->udef;
        const CHIMAT & __restrict__ CHI = pos->chi;

        for(int z=0; z<Nz; ++z)
        for(int y=0; y<Ny; ++y)
        for(int x=0; x<Nx; ++x)
        {
          // What if multiple obstacles share a block? Do not write udef onto
          // grid if CHI stored on the grid is greater than obst's CHI.
          if(c(x,y,z).s > CHI[z][y][x]) continue;
          // What if two obstacles overlap? Let's plus equal. After all here
          // we are computing divUs, maybe one obstacle has divUs 0. We will
          // need a repulsion term of the velocity at some point in the code.
          b(x,y,z).u[0] += UDEF[z][y][x][0];
          b(x,y,z).u[1] += UDEF[z][y][x][1];
          b(x,y,z).u[2] += UDEF[z][y][x][2];
        }
      }
    }
  }
}

struct KernelGradP
{
  const StencilInfo stencil{-1,-1,-1,2,2,2,false,{0}};
  SimulationData & sim;
  const std::vector<BlockInfo>& tmpVInfo = sim.tmpVInfo();
  const Real dt;
  const int Nx = VectorBlock::sizeX;
  const int Ny = VectorBlock::sizeY;
  const int Nz = VectorBlock::sizeZ;

  KernelGradP(SimulationData & s, const Real a_dt): sim(s), dt(a_dt) {}

  ~KernelGradP() {}

  void operator()(const ScalarLab & lab, const BlockInfo& info) const
  {
    VectorBlock& o = (*sim.tmpV)(info.blockID);
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
  const int Nx = VectorBlock::sizeX;
  const int Ny = VectorBlock::sizeY;
  const int Nz = VectorBlock::sizeZ;
  const std::vector<BlockInfo>& presInfo = sim.presInfo();

  pOld.resize(Nx*Ny*Nz*presInfo.size());

  //1. Compute pRHS
  {
    //pOld -> store p
    //p -> store RHS
    #pragma omp parallel for schedule(static)
    for(size_t i=0; i<presInfo.size(); i++)
    {
      const ScalarBlock& p = (*sim.pres)(i);
      VectorBlock& tmpV    = (*sim.tmpV)(i);
      for(int z=0; z<Nz; ++z)
      for(int y=0; y<Ny; ++y)
      for(int x=0; x<Nx; ++x)
      {
        pOld[i*Nx*Ny*Nz+z*Ny*Nx+y*Nx+x] = p(x,y,z).s;
        tmpV(x,y,z).u[0] = 0; 
        tmpV(x,y,z).u[1] = 0;
        tmpV(x,y,z).u[2] = 0;
      }
    }

    //place Udef on tmpV
    if(sim.obstacle_vector->nObstacles() > 0) kernelUpdateTmpV(sim);

    KernelPressureRHS K(sim,dt);
    compute<KernelPressureRHS,VectorGrid,VectorLab,VectorGrid,VectorLab,ScalarGrid>(K,*sim.vel,*sim.tmpV,true,sim.lhs);
  }

  ////2. Add div(p_old) to rhs and set initial guess phi = 0, i.e. p^{n+1}=p^{n}
  if (sim.step> sim.step_2nd_start)
  {
    compute<ScalarLab>(KernelDivPressure(sim),sim.pres,sim.tmpV);
    #pragma omp parallel for
    for(size_t i=0; i<presInfo.size(); i++)
    {
      const VectorBlock& b = (*sim.tmpV)(i);
      ScalarBlock & LHS    = (*sim.lhs )(i);
      ScalarBlock& p       = (*sim.pres)(i);
      for(int z=0; z<Nz; ++z)
      for(int y=0; y<Ny; ++y)
      for(int x=0; x<Nx; ++x)
      {
        LHS(x,y,z).s -= b(x,y,z).u[0];
        p(x,y,z).s = 0;
      }
    }
  }
  else // just set initial guess phi = 0, i.e. p^{n+1}=p^{n}
  {
    #pragma omp parallel for
    for(size_t i=0; i<presInfo.size(); i++)
    {
      ScalarBlock& p = (*sim.pres)(i);
      p.clear();
    }
  }

  //Solve the Poisson equation and use pressure to perform
  //pressure projection of the velocity field.

  //The rhs of the linear system is contained in sim.lhs
  //The initial guess is contained in sim.pres
  //Here we solve for phi := p^{n+1}-p^{n} if step < step_2nd_start
  pressureSolver->solve();

  Real avg = 0;
  Real avg1 = 0;
  #pragma omp parallel for reduction (+:avg,avg1)
  for(size_t i=0; i< presInfo.size(); i++)
  {
    ScalarBlock& P  = (*sim.pres)(i);
    const Real vv = presInfo[i].h*presInfo[i].h*presInfo[i].h;
    for(int iz=0; iz<Nz; iz++)
    for(int iy=0; iy<Ny; iy++)
    for(int ix=0; ix<Nx; ix++)
    {
      avg += P(ix,iy,iz).s * vv;
      avg1 += vv;
    }
  }
  Real quantities[2] = {avg,avg1};
  MPI_Allreduce(MPI_IN_PLACE,&quantities,2,MPI_Real,MPI_SUM,sim.comm);
  avg = quantities[0]; avg1 = quantities[1] ;
  avg = avg/avg1;
  #pragma omp parallel for
  for(size_t i=0; i< presInfo.size(); i++)
  {
    ScalarBlock & __restrict__ P  = (*sim.pres)(i);
    for(int iz=0; iz<Nz; iz++)
    for(int iy=0; iy<Ny; iy++)
    for(int ix=0; ix<Nx; ix++)
      P(ix,iy,iz).s -= avg;
  }

  const std::vector<BlockInfo>& velInfo  = sim.velInfo();

  if (sim.step > sim.step_2nd_start) //recover p^{n+1} = phi + p^{n}
  {
    #pragma omp parallel for
    for(size_t i=0; i<velInfo.size(); i++)
    {
      ScalarBlock& p          = (*sim.pres)(i);
      for (int z=0; z<Nz; ++z)
      for (int y=0; y<Ny; ++y)
      for (int x=0; x<Nx; ++x)
        p(x,y,z).s += pOld[i*Nx*Ny*Nz+z*Ny*Nx+y*Nx+x];
    }
  }

  //Compute grad(P) and put it to the vector tmpV
  compute<ScalarLab>(KernelGradP(sim,dt),sim.pres,sim.tmpV);

  //Perform the projection and set u^{n+1} = u - grad(P)*dt
  #pragma omp parallel for
  for(size_t i=0; i<velInfo.size(); i++)
  {
    const Real fac = 1.0/(velInfo[i].h*velInfo[i].h*velInfo[i].h);
    const VectorBlock& gradP = (*sim.tmpV)(i);
    VectorBlock& v = (*sim.vel)(i);
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
