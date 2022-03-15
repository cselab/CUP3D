//
//  CubismUP_3D
//  Copyright (c) 2020 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Michalis Chatzimanolakis (michaich@ethz.ch).
//

#include "PressureRHS.h"
#include "../obstacles/ObstacleVector.h"

CubismUP_3D_NAMESPACE_BEGIN
using namespace cubism;

namespace {

using CHIMAT =  Real[CUP_BLOCK_SIZEZ][CUP_BLOCK_SIZEY][CUP_BLOCK_SIZEX];
using UDEFMAT = Real[CUP_BLOCK_SIZEZ][CUP_BLOCK_SIZEY][CUP_BLOCK_SIZEX][3];

struct KernelDivPressure
{
  typedef typename FluidGridMPI::BlockType BlockType;
  const SimulationData & sim;
  const std::array<int, 3> stencil_start = {-1,-1,-1}, stencil_end = {2, 2, 2};
  const StencilInfo stencil = StencilInfo(-1,-1,-1, 2,2,2, false, {FE_P});

  KernelDivPressure(const SimulationData& s) :sim(s) {}

  template <typename Lab, typename BlockType>
  void operator()(Lab & lab, const BlockInfo& info, BlockType& o) const
  {
    BlockType & __restrict__ b  = *(BlockType*) info.ptrBlock;
    const Real fac = info.h;
    for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
    for(int iy=0; iy<FluidBlock::sizeY; ++iy)
    for(int ix=0; ix<FluidBlock::sizeX; ++ix)
      b(ix,iy,iz).tmpU = fac*(lab(ix+1,iy,iz).p+lab(ix-1,iy,iz).p
                             +lab(ix,iy+1,iz).p+lab(ix,iy-1,iz).p
                             +lab(ix,iy,iz+1).p+lab(ix,iy,iz-1).p - 6.0*lab(ix,iy,iz).p);

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
        faceXm[iy + FluidBlock::sizeY * iz].tmpU = fac *(lab(ix,iy,iz).p - lab(ix-1,iy,iz).p);
      }
    }
    if (faceXp != nullptr)
    {
      int ix = FluidBlock::sizeX-1;
      for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
      for(int iy=0; iy<FluidBlock::sizeY; ++iy)
      {
        faceXp[iy + FluidBlock::sizeY * iz].clear();
        faceXp[iy + FluidBlock::sizeY * iz].tmpU = - fac *(lab(ix+1,iy,iz).p - lab(ix,iy,iz).p);
      }
    }
    if (faceYm != nullptr)
    {
      int iy = 0;
      for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
      for(int ix=0; ix<FluidBlock::sizeX; ++ix)
      {
        faceYm[ix + FluidBlock::sizeX * iz].clear();
        faceYm[ix + FluidBlock::sizeX * iz].tmpU = fac *(lab(ix,iy,iz).p - lab(ix,iy-1,iz).p);
      }
    }
    if (faceYp != nullptr)
    {
      int iy = FluidBlock::sizeY-1;
      for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
      for(int ix=0; ix<FluidBlock::sizeX; ++ix)
      {
        faceYp[ix + FluidBlock::sizeX * iz].clear();
        faceYp[ix + FluidBlock::sizeX * iz].tmpU = - fac *(lab(ix,iy+1,iz).p - lab(ix,iy,iz).p);
      }
    }
    if (faceZm != nullptr)
    {
      int iz = 0;
      for(int iy=0; iy<FluidBlock::sizeY; ++iy)
      for(int ix=0; ix<FluidBlock::sizeX; ++ix)
      {
        faceZm[ix + FluidBlock::sizeX * iy].clear();
        faceZm[ix + FluidBlock::sizeX * iy].tmpU = fac *(lab(ix,iy,iz).p - lab(ix,iy,iz-1).p);
      }
    }
    if (faceZp != nullptr)
    {
      int iz = FluidBlock::sizeZ-1;
      for(int iy=0; iy<FluidBlock::sizeY; ++iy)
      for(int ix=0; ix<FluidBlock::sizeX; ++ix)
      {
        faceZp[ix + FluidBlock::sizeX * iy].clear();
        faceZp[ix + FluidBlock::sizeX * iy].tmpU = - fac *(lab(ix,iy,iz+1).p - lab(ix,iy,iz).p);
      }
    }
  }
};

struct KernelPressureRHS : public ObstacleVisitor
{
  typedef typename FluidGridMPI::BlockType BlockType;
  SimulationData & sim;
  const Real dt = sim.dt;
  PoissonSolverAMR * const solver = sim.pressureSolver;
  ObstacleVector * const obstacle_vector = sim.obstacle_vector;
  const int nShapes = obstacle_vector->nObstacles();
  // modified before going into accept
  const cubism::BlockInfo * info_ptr = nullptr;
  Lab * lab_ptr = nullptr;

  const std::array<int, 3> stencil_start = {-1,-1,-1}, stencil_end = {2, 2, 2};
  const StencilInfo stencil = StencilInfo(-1,-1,-1, 2,2,2, false, {FE_U, FE_V, FE_W, FE_TMPU, FE_TMPV, FE_TMPW});

  KernelPressureRHS(SimulationData& s) :sim(s) {}

  template <typename Lab, typename BlockType>
  void operator()(Lab & lab, const BlockInfo& info, BlockType& o)
  {
    const Real h = info.h, fac = 0.5*h*h/dt;
    BlockType & __restrict__ b  = *(BlockType*) info.ptrBlock;

    for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
    for(int iy=0; iy<FluidBlock::sizeY; ++iy)
    for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
      const FluidElement &LW = lab(ix-1,iy,  iz  ), &LE = lab(ix+1,iy,  iz  );
      const FluidElement &LS = lab(ix,  iy-1,iz  ), &LN = lab(ix,  iy+1,iz  );
      const FluidElement &LF = lab(ix,  iy,  iz-1), &LB = lab(ix,  iy,  iz+1);
      b(ix,iy,iz).p = fac*(LE.u-LW.u + LN.v-LS.v + LB.w-LF.w);
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
        faceXm[iy + FluidBlock::sizeY * iz].p = fac *(lab(ix-1,iy,iz).u + lab(ix,iy,iz).u);
      }
    }
    if (faceXp != nullptr)
    {
      int ix = FluidBlock::sizeX-1;
      for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
      for(int iy=0; iy<FluidBlock::sizeY; ++iy)
      {
        faceXp[iy + FluidBlock::sizeY * iz].clear();
        faceXp[iy + FluidBlock::sizeY * iz].p = - fac *(lab(ix+1,iy,iz).u + lab(ix,iy,iz).u);
      }
    }
    if (faceYm != nullptr)
    {
      int iy = 0;
      for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
      for(int ix=0; ix<FluidBlock::sizeX; ++ix)
      {
        faceYm[ix + FluidBlock::sizeX * iz].clear();
        faceYm[ix + FluidBlock::sizeX * iz].p = fac *(lab(ix,iy-1,iz).v + lab(ix,iy,iz).v);
      }
    }
    if (faceYp != nullptr)
    {
      int iy = FluidBlock::sizeY-1;
      for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
      for(int ix=0; ix<FluidBlock::sizeX; ++ix)
      {
        faceYp[ix + FluidBlock::sizeX * iz].clear();
        faceYp[ix + FluidBlock::sizeX * iz].p = - fac *(lab(ix,iy+1,iz).v + lab(ix,iy,iz).v);
      }
    }
    if (faceZm != nullptr)
    {
      int iz = 0;
      for(int iy=0; iy<FluidBlock::sizeY; ++iy)
      for(int ix=0; ix<FluidBlock::sizeX; ++ix)
      {
        faceZm[ix + FluidBlock::sizeX * iy].clear();
        faceZm[ix + FluidBlock::sizeX * iy].p = fac *(lab(ix,iy,iz-1).w + lab(ix,iy,iz).w);
      }
    }
    if (faceZp != nullptr)
    {
      int iz = FluidBlock::sizeZ-1;
      for(int iy=0; iy<FluidBlock::sizeY; ++iy)
      for(int ix=0; ix<FluidBlock::sizeX; ++ix)
      {
        faceZp[ix + FluidBlock::sizeX * iy].clear();
        faceZp[ix + FluidBlock::sizeX * iy].p = - fac *(lab(ix,iy,iz+1).w + lab(ix,iy,iz).w);
      }
    }
    if(nShapes == 0) return; // no need to account for obstacles

    // first store the lab and info, then do visitor
    assert(info_ptr == nullptr && lab_ptr == nullptr);
    info_ptr =  & info; lab_ptr =   & lab;
    ObstacleVisitor* const base = static_cast<ObstacleVisitor*> (this);
    assert( base not_eq nullptr );
    obstacle_vector->Accept( base );
    info_ptr = nullptr; lab_ptr = nullptr;
  }

  void visit(Obstacle* const obstacle)
  {
    assert(info_ptr not_eq nullptr && lab_ptr not_eq nullptr);
    const BlockInfo& info = * info_ptr;
    Lab& lab = * lab_ptr;
    const auto& obstblocks = obstacle->getObstacleBlocks();
    if (obstblocks[info.blockID] == nullptr) return;

    const CHIMAT & __restrict__ CHI = obstblocks[info.blockID]->chi;

    const Real h = info.h, fac = 0.5*h*h/dt;

    BlockType & __restrict__ b  = *(BlockType*) info.ptrBlock;
    for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
    for(int iy=0; iy<FluidBlock::sizeY; ++iy)
    for(int ix=0; ix<FluidBlock::sizeX; ++ix)
    {
      if (lab(ix,iy,iz).chi > CHI[iz][iy][ix]) continue;
      const FluidElement &LW = lab(ix-1,iy,  iz  ), &LE = lab(ix+1,iy,  iz  );
      const FluidElement &LS = lab(ix,  iy-1,iz  ), &LN = lab(ix,  iy+1,iz  );
      const FluidElement &LF = lab(ix,  iy,  iz-1), &LB = lab(ix,  iy,  iz+1);
      const Real divUs = LE.tmpU-LW.tmpU + LN.tmpV-LS.tmpV + LB.tmpW-LF.tmpW;
      const Real srcBulk = - CHI[iz][iy][ix] * fac * divUs;
      b(ix,iy,iz).p += srcBulk;
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
        faceXm[iy + FluidBlock::sizeY * iz].p += - CHI[iz][iy][ix]* fac *(lab(ix-1,iy,iz).tmpU + lab(ix,iy,iz).tmpU);
      }
    }
    if (faceXp != nullptr)
    {
      int ix = FluidBlock::sizeX-1;
      for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
      for(int iy=0; iy<FluidBlock::sizeY; ++iy)
      {
        faceXp[iy + FluidBlock::sizeY * iz].p +=   CHI[iz][iy][ix] * fac *(lab(ix+1,iy,iz).tmpU + lab(ix,iy,iz).tmpU);
      }
    }
    if (faceYm != nullptr)
    {
      int iy = 0;
      for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
      for(int ix=0; ix<FluidBlock::sizeX; ++ix)
      {
        faceYm[ix + FluidBlock::sizeX * iz].p += - CHI[iz][iy][ix] * fac *(lab(ix,iy-1,iz).tmpV + lab(ix,iy,iz).tmpV);
      }
    }
    if (faceYp != nullptr)
    {
      int iy = FluidBlock::sizeY-1;
      for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
      for(int ix=0; ix<FluidBlock::sizeX; ++ix)
      {
        faceYp[ix + FluidBlock::sizeX * iz].p +=  CHI[iz][iy][ix]* fac *(lab(ix,iy+1,iz).tmpV + lab(ix,iy,iz).tmpV);
      }
    }
    if (faceZm != nullptr)
    {
      int iz = 0;
      for(int iy=0; iy<FluidBlock::sizeY; ++iy)
      for(int ix=0; ix<FluidBlock::sizeX; ++ix)
      {
        faceZm[ix + FluidBlock::sizeX * iy].p += - CHI[iz][iy][ix] *fac *(lab(ix,iy,iz-1).tmpW + lab(ix,iy,iz).tmpW);
      }
    }
    if (faceZp != nullptr)
    {
      int iz = FluidBlock::sizeZ-1;
      for(int iy=0; iy<FluidBlock::sizeY; ++iy)
      for(int ix=0; ix<FluidBlock::sizeX; ++ix)
      {
        faceZp[ix + FluidBlock::sizeX * iy].p +=  CHI[iz][iy][ix]*fac *(lab(ix,iy,iz+1).tmpW + lab(ix,iy,iz).tmpW);
      }
    }
  }
};

struct PressureRHSObstacleVisitor : public ObstacleVisitor
{
  typedef typename FluidGridMPI::BlockType BlockType;
  FluidGridMPI * const grid;
  const std::vector<cubism::BlockInfo>& vInfo = grid->getBlocksInfo();

  PressureRHSObstacleVisitor(FluidGridMPI*g) : grid(g) { }

  void visit(Obstacle* const obstacle)
  {
    #pragma omp parallel
    {
      const auto& obstblocks = obstacle->getObstacleBlocks();
      #pragma omp for schedule(dynamic, 1)
      for (size_t i = 0; i < vInfo.size(); ++i)
      {
        const cubism::BlockInfo& info = vInfo[i];
        const auto pos = obstblocks[info.blockID];
        if(pos == nullptr) continue;

        FluidBlock& b = *(FluidBlock*)info.ptrBlock;
        const UDEFMAT & __restrict__ UDEF = pos->udef;
        const CHIMAT & __restrict__ CHI = pos->chi;

        for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
        for(int iy=0; iy<FluidBlock::sizeY; ++iy)
        for(int ix=0; ix<FluidBlock::sizeX; ++ix)
        {
          // What if multiple obstacles share a block? Do not write udef onto
          // grid if CHI stored on the grid is greater than obst's CHI.
          if(b(ix,iy,iz).chi > CHI[iz][iy][ix]) continue;
          // What if two obstacles overlap? Let's plus equal. After all here
          // we are computing divUs, maybe one obstacle has divUs 0. We will
          // need a repulsion term of the velocity at some point in the code.
          b(ix,iy,iz).tmpU += UDEF[iz][iy][ix][0];
          b(ix,iy,iz).tmpV += UDEF[iz][iy][ix][1];
          b(ix,iy,iz).tmpW += UDEF[iz][iy][ix][2];
        }
      }
    }
  }
};

}

PressureRHS::PressureRHS(SimulationData & s) : Operator(s) {}

void PressureRHS::operator()(const Real dt)
{
  const std::vector<cubism::BlockInfo>& vInfo = grid->getBlocksInfo();
  const std::vector<cubism::BlockInfo>& vInfoPoisson = gridPoisson->getBlocksInfo();

  //1. Compute pRHS
  {
    //dataOld[3] -> store p
    //p -> store RHS
    #pragma omp parallel for schedule(static)
    for(size_t i=0; i<vInfo.size(); i++)
    {
      FluidBlock& b = *(FluidBlock*)vInfo[i].ptrBlock;
      for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
      for(int iy=0; iy<FluidBlock::sizeY; ++iy)
      for(int ix=0; ix<FluidBlock::sizeX; ++ix)
      {
	b.dataOld[iz][iy][ix][3] = b(ix,iy,iz).p;//store old pressure where it should be stored
        b(ix,iy,iz).tmpU = 0; 
	b(ix,iy,iz).tmpV = 0;
	b(ix,iy,iz).tmpW = 0;
      }
    }

    //place Udef on tmpU,tmpV,tmpW
    const size_t nShapes = sim.obstacle_vector->nObstacles();
    if(nShapes > 0)
    {
      ObstacleVisitor* visitor = new PressureRHSObstacleVisitor(grid);
      sim.obstacle_vector->Accept(visitor);
      delete visitor;
    }

    KernelPressureRHS K(sim);
    compute<KernelPressureRHS> (K,true);

    #pragma omp parallel for schedule(static)
    for(size_t i=0; i<vInfo.size(); i++)
    {
      FluidBlock& b = *(FluidBlock*)vInfo[i].ptrBlock;
      FluidBlockPoisson& bPoisson = *(FluidBlockPoisson*)vInfoPoisson[i].ptrBlock;
      for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
      for(int iy=0; iy<FluidBlock::sizeY; ++iy)
      for(int ix=0; ix<FluidBlock::sizeX; ++ix)
      {
        bPoisson(ix,iy,iz).lhs = b(ix,iy,iz).p;
        b(ix,iy,iz).p = b.dataOld[iz][iy][ix][3];
      }
    }
  }

  ////2. Add div(p_old) to rhs
  if (sim.step> sim.step_2nd_start)
  {
    const KernelDivPressure K(sim);
    compute(K,true);
    #pragma omp parallel for
    for(size_t i=0; i<vInfo.size(); i++)
    {
      FluidBlock& b = *(FluidBlock*)vInfo[i].ptrBlock;
      FluidBlockPoisson& bPoisson = *(FluidBlockPoisson*)vInfoPoisson[i].ptrBlock;
      for (int iz=0; iz<FluidBlock::sizeZ; ++iz)
      for (int iy=0; iy<FluidBlock::sizeY; ++iy)
      for (int ix=0; ix<FluidBlock::sizeX; ++ix)
        bPoisson(ix,iy,iz).lhs -= b(ix,iy,iz).tmpU;
    }
  }

  //initial guess phi = 0, i.e. p^{n+1}=p^{n}
  #pragma omp parallel for
  for(size_t i=0; i<vInfo.size(); i++)
  {
    FluidBlock& b = *(FluidBlock*)vInfo[i].ptrBlock;
    FluidBlockPoisson& bPoisson = *(FluidBlockPoisson*)vInfoPoisson[i].ptrBlock;
    for (int iz=0; iz<FluidBlock::sizeZ; ++iz)
    for (int iy=0; iy<FluidBlock::sizeY; ++iy)
    for (int ix=0; ix<FluidBlock::sizeX; ++ix)
    {
      bPoisson(ix,iy,iz).s = 0.0;//pass initial guess to solver
      b(ix,iy,iz).p = 0.0;
    }
  }
  check("PressureRHS");
}

CubismUP_3D_NAMESPACE_END
