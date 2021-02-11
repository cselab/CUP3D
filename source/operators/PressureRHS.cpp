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

static constexpr Real EPS = std::numeric_limits<Real>::epsilon();
using CHIMAT = Real[CUP_BLOCK_SIZE][CUP_BLOCK_SIZE][CUP_BLOCK_SIZE];
using UDEFMAT = Real[CUP_BLOCK_SIZE][CUP_BLOCK_SIZE][CUP_BLOCK_SIZE][3];

struct KernelPressureRHS
{
  typedef typename FluidGridMPI::BlockType BlockType;
  SimulationData & sim;
  const Real dt = sim.dt;
  PoissonSolverAMR * const solver = sim.pressureSolver;
  ObstacleVector * const obstacle_vector = sim.obstacle_vector;
  const int nShapes = obstacle_vector->nObstacles();
  // non-const non thread safe:
  std::vector<Real> sumRHS = std::vector<Real>(nShapes, 0);
  std::vector<Real> posRHS = std::vector<Real>(nShapes, 0);
  std::vector<Real> negRHS = std::vector<Real>(nShapes, 0);
  // modified before going into accept
  const cubism::BlockInfo * info_ptr = nullptr;
  Lab * lab_ptr = nullptr;

  const std::array<int, 3> stencil_start = {-1,-1,-1}, stencil_end = {2, 2, 2};
  const StencilInfo stencil = StencilInfo(-1,-1,-1, 2,2,2, false, {FE_U, FE_V, FE_W, FE_TMPU, FE_TMPV, FE_TMPW});

  KernelPressureRHS(SimulationData& s) :sim(s) {}

  template <typename Lab, typename BlockType>
  void operator()(Lab & lab, const BlockInfo& info, BlockType& o)
  {
    const Real h = info.h_gridpoint, fac = 0.5*h*h/dt;
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
  }
};

}

PressureRHS::PressureRHS(SimulationData & s) : Operator(s) {}

void PressureRHS::operator()(const double dt)
{
  sim.startProfiler("PresRHS Kernel");

  sim.pressureSolver->reset();

  const std::vector<cubism::BlockInfo>& vInfo = grid->getBlocksInfo();
  const int nthreads = omp_get_max_threads();

  //tmp -> store p
  //p -> store RHS
  #pragma omp parallel for schedule(static)
  for(size_t i=0; i<vInfo.size(); i++) 
  {
    FluidBlock& b = *(FluidBlock*)vInfo[i].ptrBlock;
    for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
    for(int iy=0; iy<FluidBlock::sizeY; ++iy)
    for(int ix=0; ix<FluidBlock::sizeX; ++ix)
      std::swap(b(ix,iy,iz).p, b.tmp[iz][iy][ix]);
  }

  std::vector< KernelPressureRHS *> K(nthreads, nullptr);
  for(int i=0;i<nthreads;++i) K[i] = new KernelPressureRHS (sim);
  compute< KernelPressureRHS  >(K,true);//true: apply FluxCorrection
  for(int i=0; i<nthreads; i++) delete K[i];

  //tmp -> store RHS
  //(p and tmp swap)
  #pragma omp parallel for schedule(static)
  for(size_t i=0; i<vInfo.size(); i++) 
  {
    FluidBlock& b = *(FluidBlock*)vInfo[i].ptrBlock;
    for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
    for(int iy=0; iy<FluidBlock::sizeY; ++iy)
    for(int ix=0; ix<FluidBlock::sizeX; ++ix)
      std::swap(b(ix,iy,iz).p, b.tmp[iz][iy][ix]);
  }

  check("PressureRHS");
  sim.stopProfiler();
}

CubismUP_3D_NAMESPACE_END
