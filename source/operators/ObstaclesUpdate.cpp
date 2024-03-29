//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#include "ObstaclesUpdate.h"
#include "../Obstacles/ObstacleVector.h"
#include "../Utils/MatArrayMath.h"

// define this to update obstacles with old (mrag-like) approach of integrating
// momenta contained in chi before the penalization step:

CubismUP_3D_NAMESPACE_BEGIN
using namespace cubism;

namespace {

using CHIMAT =  Real[CUP_BLOCK_SIZEZ][CUP_BLOCK_SIZEY][CUP_BLOCK_SIZEX];
using UDEFMAT = Real[CUP_BLOCK_SIZEZ][CUP_BLOCK_SIZEY][CUP_BLOCK_SIZEX][3];

template<bool implicitPenalization>
struct KernelIntegrateFluidMomenta
{
  const Real lambda, dt;
  ObstacleVector * const obstacle_vector;

  Real dvol(const BlockInfo&I, const int x, const int y, const int z) const {
    return I.h * I.h * I.h;
  }

  KernelIntegrateFluidMomenta(Real _dt, Real _lambda, ObstacleVector* ov)
    : lambda(_lambda), dt(_dt), obstacle_vector(ov) {}

  void operator()(const cubism::BlockInfo& info) const
  {
    for (const auto &obstacle : obstacle_vector->getObstacleVector())
      visit(info, obstacle.get());
  }

  void visit(const BlockInfo& info, Obstacle* const op) const
  {
    const std::vector<ObstacleBlock*>& obstblocks = op->getObstacleBlocks();
    ObstacleBlock*const o = obstblocks[info.blockID];
    if (o == nullptr) return;

    const std::array<Real,3> CM = op->getCenterOfMass();
    const VectorBlock &b = *(VectorBlock *)info.ptrBlock;
    const CHIMAT & __restrict__ CHI = o->chi;
    Real &VV = o->V;
    Real &FX = o->FX, &FY = o->FY, &FZ = o->FZ;
    Real &TX = o->TX, &TY = o->TY, &TZ = o->TZ;
    VV = 0; FX = 0; FY = 0; FZ = 0; TX = 0; TY = 0; TZ = 0;
    Real &J0 = o->J0, &J1 = o->J1, &J2 = o->J2;
    Real &J3 = o->J3, &J4 = o->J4, &J5 = o->J5;
    J0 = 0; J1 = 0; J2 = 0; J3 = 0; J4 = 0; J5 = 0;

    const UDEFMAT & __restrict__ UDEF = o->udef;
    const Real lambdt = lambda*dt;
    if(implicitPenalization)
    {
      o->GfX = 0;
      o->GpX = 0; o->GpY = 0; o->GpZ = 0;
      o->Gj0 = 0; o->Gj1 = 0; o->Gj2 = 0;
      o->Gj3 = 0; o->Gj4 = 0; o->Gj5 = 0;
      o->GuX = 0; o->GuY = 0; o->GuZ = 0;
      o->GaX = 0; o->GaY = 0; o->GaZ = 0;
    }

    for(int iz=0; iz<VectorBlock::sizeZ; ++iz)
    for(int iy=0; iy<VectorBlock::sizeY; ++iy)
    for(int ix=0; ix<VectorBlock::sizeX; ++ix)
    {
      if (CHI[iz][iy][ix] <= 0) continue;
      Real p[3]; info.pos(p, ix, iy, iz);
      const Real dv = dvol(info, ix, iy, iz), X = CHI[iz][iy][ix];
      p[0] -= CM[0]; p[1] -= CM[1]; p[2] -= CM[2];

      VV += X * dv;
      J0 += X * dv * ( p[1]*p[1] + p[2]*p[2] );
      J1 += X * dv * ( p[0]*p[0] + p[2]*p[2] );
      J2 += X * dv * ( p[0]*p[0] + p[1]*p[1] );
      J3 -= X * dv * p[0]*p[1];
      J4 -= X * dv * p[0]*p[2];
      J5 -= X * dv * p[1]*p[2];

      FX += X * dv * b(ix,iy,iz).u[0];
      FY += X * dv * b(ix,iy,iz).u[1];
      FZ += X * dv * b(ix,iy,iz).u[2];
      TX += X * dv * ( p[1] * b(ix,iy,iz).u[2] - p[2] * b(ix,iy,iz).u[1] );
      TY += X * dv * ( p[2] * b(ix,iy,iz).u[0] - p[0] * b(ix,iy,iz).u[2] );
      TZ += X * dv * ( p[0] * b(ix,iy,iz).u[1] - p[1] * b(ix,iy,iz).u[0] );

      if(implicitPenalization)
      {
        const Real X1 = CHI[iz][iy][ix]>0.5?1.0:0.0;
        const Real penalFac = dv * lambdt * X1 / ( 1 +  X1 * lambdt );
        o->GfX += penalFac;
        o->GpX += penalFac * p[0];
        o->GpY += penalFac * p[1];
        o->GpZ += penalFac * p[2];
        o->Gj0 += penalFac * ( p[1]*p[1] + p[2]*p[2] );
        o->Gj1 += penalFac * ( p[0]*p[0] + p[2]*p[2] );
        o->Gj2 += penalFac * ( p[0]*p[0] + p[1]*p[1] );
        o->Gj3 -= penalFac * p[0]*p[1];
        o->Gj4 -= penalFac * p[0]*p[2];
        o->Gj5 -= penalFac * p[1]*p[2];
        const Real DiffU[3] = {
          b(ix,iy,iz).u[0] - UDEF[iz][iy][ix][0],
          b(ix,iy,iz).u[1] - UDEF[iz][iy][ix][1],
          b(ix,iy,iz).u[2] - UDEF[iz][iy][ix][2]
        };
        o->GuX += penalFac * DiffU[0];
        o->GuY += penalFac * DiffU[1];
        o->GuZ += penalFac * DiffU[2];
        o->GaX += penalFac * ( p[1] * DiffU[2] - p[2] * DiffU[1] );
        o->GaY += penalFac * ( p[2] * DiffU[0] - p[0] * DiffU[2] );
        o->GaZ += penalFac * ( p[0] * DiffU[1] - p[1] * DiffU[0] );
      }
    }
  }
};

}  // Anonymous namespace.

template<bool implicitPenalization>
static void kernelFinalizeObstacleVel(SimulationData& sim, const Real dt)
{
  // TODO: Refactor to use only one omp parallel and one MPI_Allreduce.
  for (const auto &obst : sim.obstacle_vector->getObstacleVector()) {
    static constexpr int nQoI = 29;
    Real M[nQoI] = { 0 };
    const auto& oBlock = obst->getObstacleBlocks();
    #pragma omp parallel for schedule(static,1) reduction(+ : M[:nQoI])
    for (size_t i=0; i<oBlock.size(); i++) {
      if(oBlock[i] == nullptr) continue;
      int k = 0;
      M[k++] += oBlock[i]->V ;
      M[k++] += oBlock[i]->FX; M[k++] += oBlock[i]->FY; M[k++] += oBlock[i]->FZ;
      M[k++] += oBlock[i]->TX; M[k++] += oBlock[i]->TY; M[k++] += oBlock[i]->TZ;
      M[k++] += oBlock[i]->J0; M[k++] += oBlock[i]->J1; M[k++] += oBlock[i]->J2;
      M[k++] += oBlock[i]->J3; M[k++] += oBlock[i]->J4; M[k++] += oBlock[i]->J5;
      if(implicitPenalization) {
      M[k++] +=oBlock[i]->GfX;
      M[k++] +=oBlock[i]->GpX; M[k++] +=oBlock[i]->GpY; M[k++] +=oBlock[i]->GpZ;
      M[k++] +=oBlock[i]->Gj0; M[k++] +=oBlock[i]->Gj1; M[k++] +=oBlock[i]->Gj2;
      M[k++] +=oBlock[i]->Gj3; M[k++] +=oBlock[i]->Gj4; M[k++] +=oBlock[i]->Gj5;
      M[k++] +=oBlock[i]->GuX; M[k++] +=oBlock[i]->GuY; M[k++] +=oBlock[i]->GuZ;
      M[k++] +=oBlock[i]->GaX; M[k++] +=oBlock[i]->GaY; M[k++] +=oBlock[i]->GaZ;
      assert(k==29);
      } else  assert(k==13);
    }
    const auto comm = sim.comm;
    MPI_Allreduce(MPI_IN_PLACE, M, nQoI, MPI_Real, MPI_SUM, comm);

    #ifndef NDEBUG
      const Real J_magnitude = obst->J[0] + obst->J[1] + obst->J[2];
      static constexpr Real EPS = std::numeric_limits<Real>::epsilon();
    #endif
    assert(std::fabs(obst->mass - M[ 0]) < 10 * EPS * obst->mass);
    assert(std::fabs(obst->J[0] - M[ 7]) < 10 * EPS * J_magnitude);
    assert(std::fabs(obst->J[1] - M[ 8]) < 10 * EPS * J_magnitude);
    assert(std::fabs(obst->J[2] - M[ 9]) < 10 * EPS * J_magnitude);
    assert(std::fabs(obst->J[3] - M[10]) < 10 * EPS * J_magnitude);
    assert(std::fabs(obst->J[4] - M[11]) < 10 * EPS * J_magnitude);
    assert(std::fabs(obst->J[5] - M[12]) < 10 * EPS * J_magnitude);
    assert(M[0] > EPS);

    if(implicitPenalization) {
      obst->penalM    = M[13];
      obst->penalCM   = { M[14], M[15], M[16] };
      obst->penalJ    = { M[17], M[18], M[19], M[20], M[21], M[22] };
      obst->penalLmom = { M[23], M[24], M[25] };
      obst->penalAmom = { M[26], M[27], M[28] };
    } else {
      obst->penalM    = M[0];
      obst->penalCM   = { 0, 0, 0 };
      obst->penalJ    = { M[ 7], M[ 8], M[ 9], M[10], M[11], M[12] };
      obst->penalLmom = { M[1], M[2], M[3] };
      obst->penalAmom = { M[4], M[5], M[6] };
    }

    obst->computeVelocities();
  }
}

void UpdateObstacles::operator()(const Real dt)
{
  if(sim.obstacle_vector->nObstacles() == 0) return;

  { // integrate momenta by looping over grid
    std::vector<cubism::BlockInfo>& velInfo = sim.velInfo();
    #pragma omp parallel
    {
      //if(0) {
      if(sim.bImplicitPenalization) {
        KernelIntegrateFluidMomenta<1> K(dt, sim.lambda, sim.obstacle_vector);
        #pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < velInfo.size(); ++i) K(velInfo[i]);
      } else {
        KernelIntegrateFluidMomenta<0> K(dt, sim.lambda, sim.obstacle_vector);
        #pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < velInfo.size(); ++i) K(velInfo[i]);
      }
    }
  }

  //if(0) {
  if(sim.bImplicitPenalization) {
    kernelFinalizeObstacleVel<1>(sim, dt);
  } else {
    kernelFinalizeObstacleVel<0>(sim, dt);
  }
}

CubismUP_3D_NAMESPACE_END
