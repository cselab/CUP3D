//
//  CubismUP_3D
//  Copyright (c) 2020 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Michalis Chatzimanolakis (michaich@ethz.ch).
//

#include "PressureProjection.h"

CubismUP_3D_NAMESPACE_BEGIN
using namespace cubism;

namespace {

class KernelGradP
{
  const Real dt;
 public:
  const std::array<int, 3> stencil_start = {-1,-1,-1}, stencil_end = {2, 2, 2};
  const StencilInfo stencil{-1,-1,-1, 2,2,2, false, {FE_P} };

  KernelGradP(double _dt): dt(_dt) {}

  ~KernelGradP() {}

  template <typename Lab, typename BlockType>
  void operator()(Lab & lab, const BlockInfo& info, BlockType& o) const
  {
    const Real fac = -0.5*dt*info.h*info.h;
    for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
    for(int iy=0; iy<FluidBlock::sizeY; ++iy)
    for(int ix=0; ix<FluidBlock::sizeX; ++ix)
    {
      o(ix,iy,iz).tmpU = fac*(lab(ix+1,iy,iz).p-lab(ix-1,iy,iz).p);
      o(ix,iy,iz).tmpV = fac*(lab(ix,iy+1,iz).p-lab(ix,iy-1,iz).p);
      o(ix,iy,iz).tmpW = fac*(lab(ix,iy,iz+1).p-lab(ix,iy,iz-1).p);
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
        faceXm[iy + FluidBlock::sizeY * iz].tmpU = fac *(lab(ix-1,iy,iz).p + lab(ix,iy,iz).p);
      }
    }
    if (faceXp != nullptr)
    {
      int ix = FluidBlock::sizeX-1;
      for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
      for(int iy=0; iy<FluidBlock::sizeY; ++iy)
      {
        faceXp[iy + FluidBlock::sizeY * iz].clear();
        faceXp[iy + FluidBlock::sizeY * iz].tmpU = - fac *(lab(ix+1,iy,iz).p + lab(ix,iy,iz).p);
      }
    }
    if (faceYm != nullptr)
    {
      int iy = 0;
      for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
      for(int ix=0; ix<FluidBlock::sizeX; ++ix)
      {
        faceYm[ix + FluidBlock::sizeX * iz].clear();
        faceYm[ix + FluidBlock::sizeX * iz].tmpV = fac *(lab(ix,iy-1,iz).p + lab(ix,iy,iz).p);
      }
    }
    if (faceYp != nullptr)
    {
      int iy = FluidBlock::sizeY-1;
      for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
      for(int ix=0; ix<FluidBlock::sizeX; ++ix)
      {
        faceYp[ix + FluidBlock::sizeX * iz].clear();
        faceYp[ix + FluidBlock::sizeX * iz].tmpV = - fac *(lab(ix,iy+1,iz).p + lab(ix,iy,iz).p);
      }
    }
    if (faceZm != nullptr)
    {
      int iz = 0;
      for(int iy=0; iy<FluidBlock::sizeY; ++iy)
      for(int ix=0; ix<FluidBlock::sizeX; ++ix)
      {
        faceZm[ix + FluidBlock::sizeX * iy].clear();
        faceZm[ix + FluidBlock::sizeX * iy].tmpW = fac *(lab(ix,iy,iz-1).p + lab(ix,iy,iz).p);
      }
    }
    if (faceZp != nullptr)
    {
      int iz = FluidBlock::sizeZ-1;
      for(int iy=0; iy<FluidBlock::sizeY; ++iy)
      for(int ix=0; ix<FluidBlock::sizeX; ++ix)
      {
        faceZp[ix + FluidBlock::sizeX * iy].clear();
        faceZp[ix + FluidBlock::sizeX * iy].tmpW = - fac *(lab(ix,iy,iz+1).p + lab(ix,iy,iz).p);
      }
    }
  }
};

}

PressureProjection::PressureProjection(SimulationData & s) : Operator(s)
{
  pressureSolver = new PoissonSolverAMR(sim);
  sim.pressureSolver = pressureSolver;
}

void PressureProjection::operator()(const double dt)
{
  //The initial guess is contained in vInfoPoisson -> s
  //The rhs is contained in vInfoPoisson -> lhs
  const std::vector<cubism::BlockInfo>& vInfo = sim.vInfo();

  // solve for phi := p^{n+1}-p^{n}
  pressureSolver->solve();//will return p=phi
  if (sim.step > sim.step_2nd_start) //recover p^{n+1} from phi
  {
    #pragma omp parallel for
    for(size_t i=0; i<vInfo.size(); i++)
    {
      FluidBlock& b = *(FluidBlock*)vInfo[i].ptrBlock;
      for (int iz=0; iz<FluidBlock::sizeZ; ++iz)
      for (int iy=0; iy<FluidBlock::sizeY; ++iy)
      for (int ix=0; ix<FluidBlock::sizeX; ++ix)
        b(ix,iy,iz).p += b.dataOld[iz][iy][ix][3];
    }
  }

  //pressure correction dudt* = - grad P / rho
  const KernelGradP K(dt);
  compute<KernelGradP>(K,true);

  #pragma omp parallel for
  for(size_t i=0; i<vInfo.size(); i++)
  {
    const double fac = 1.0/(vInfo[i].h*vInfo[i].h*vInfo[i].h);
    FluidBlock& b = *(FluidBlock*)vInfo[i].ptrBlock;
    for (int iz=0; iz<FluidBlock::sizeZ; ++iz)
    for (int iy=0; iy<FluidBlock::sizeY; ++iy)
    for (int ix=0; ix<FluidBlock::sizeX; ++ix)
    {
      b(ix,iy,iz).u += fac*b(ix,iy,iz).tmpU;
      b(ix,iy,iz).v += fac*b(ix,iy,iz).tmpV;
      b(ix,iy,iz).w += fac*b(ix,iy,iz).tmpW;
    }
  }

  check("PressureProjection");
}

CubismUP_3D_NAMESPACE_END
