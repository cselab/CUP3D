//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//

#include "ObstaclesCreate.h"

CubismUP_3D_NAMESPACE_BEGIN
using namespace cubism;

namespace {

using CHIMAT =  Real[CUP_BLOCK_SIZEZ][CUP_BLOCK_SIZEY][CUP_BLOCK_SIZEX];
using UDEFMAT = Real[CUP_BLOCK_SIZEZ][CUP_BLOCK_SIZEY][CUP_BLOCK_SIZEX][3];
static constexpr Real EPS = std::numeric_limits<Real>::epsilon();

struct KernelCharacteristicFunction
{
  using v_v_ob = std::vector<std::vector<ObstacleBlock*>*>;
  const v_v_ob & vec_obstacleBlocks;
  const int Nx = VectorBlock::sizeX;
  const int Ny = VectorBlock::sizeY;
  const int Nz = VectorBlock::sizeZ;

  KernelCharacteristicFunction(const v_v_ob& v) : vec_obstacleBlocks(v) {}

  void operate(const BlockInfo & info, ScalarBlock & b) const
  {
    const Real h = info.h, inv2h = .5/h, fac1 = .5*h*h, vol = h*h*h;
    const int gp = 1;

    for (size_t obst_id = 0; obst_id<vec_obstacleBlocks.size(); obst_id++)
    {
      const auto& obstacleBlocks = * vec_obstacleBlocks[obst_id];
      ObstacleBlock* const o = obstacleBlocks[info.blockID];
      if(o == nullptr) continue;
      CHIMAT & __restrict__ CHI = o->chi;
      o->CoM_x = 0; o->CoM_y = 0; o->CoM_z = 0; o->mass  = 0;
      const auto & SDFLAB = o->sdfLab;
      //////////////////////////
      // FDMH_1 computation to approximate Heaviside function H(SDF(x,y,z))
      // Reference: John D.Towers, "Finite difference methods for approximating Heaviside functions", eq.(14)
      //////////////////////////
      for(int z=0; z<Nz; ++z)
      for(int y=0; y<Ny; ++y)
      for(int x=0; x<Nx; ++x)
      {
        #if 1
        // here I read fist from SDF to deal with obstacles sharing block
        if (SDFLAB[z+1][y+1][x+1] > +gp*h || SDFLAB[z+1][y+1][x+1] < -gp*h)
        {
          CHI[z][y][x] = SDFLAB[z+1][y+1][x+1] > 0 ? 1 : 0;
        }
        else
        {
          const Real distPx = SDFLAB[z+1][y+1][x+1+1];
          const Real distMx = SDFLAB[z+1][y+1][x+1-1];
          const Real distPy = SDFLAB[z+1][y+1+1][x+1];
          const Real distMy = SDFLAB[z+1][y+1-1][x+1];
          const Real distPz = SDFLAB[z+1+1][y+1][x+1];
          const Real distMz = SDFLAB[z+1-1][y+1][x+1];
          // gradU
          const Real gradUX = inv2h*(distPx - distMx);
          const Real gradUY = inv2h*(distPy - distMy);
          const Real gradUZ = inv2h*(distPz - distMz);
          const Real gradUSq = gradUX*gradUX+gradUY*gradUY+gradUZ*gradUZ + EPS;
          const Real IplusX = std::max((Real)0.0,distPx);
          const Real IminuX = std::max((Real)0.0,distMx);
          const Real IplusY = std::max((Real)0.0,distPy);
          const Real IminuY = std::max((Real)0.0,distMy);
          const Real IplusZ = std::max((Real)0.0,distPz);
          const Real IminuZ = std::max((Real)0.0,distMz);
          // gradI: first primitive of H(x): I(x) = int_0^x H(y) dy
          const Real gradIX = inv2h*(IplusX - IminuX);
          const Real gradIY = inv2h*(IplusY - IminuY);
          const Real gradIZ = inv2h*(IplusZ - IminuZ);
          const Real numH = gradIX*gradUX + gradIY*gradUY + gradIZ*gradUZ;
          CHI[z][y][x] = numH/gradUSq;
        }
        #else
          CHI[z][y][x] = SDFLAB[z+1][y+1][x+1] > 0 ? 1 : 0;
        #endif
        Real p[3]; info.pos(p,x,y,z);
        b(x,y,z).s = std::max(CHI[z][y][x], b(x,y,z).s);
        o->CoM_x += CHI[z][y][x] * vol * p[0];
        o->CoM_y += CHI[z][y][x] * vol * p[1];
        o->CoM_z += CHI[z][y][x] * vol * p[2];
        o->mass  += CHI[z][y][x] * vol;
      }

      for(int z=0; z<Nz; ++z)
      for(int y=0; y<Ny; ++y)
      for(int x=0; x<Nx; ++x)
      {
          const Real distPx = SDFLAB[z+1][y+1][x+1+1];
          const Real distMx = SDFLAB[z+1][y+1][x+1-1];
          const Real distPy = SDFLAB[z+1][y+1+1][x+1];
          const Real distMy = SDFLAB[z+1][y+1-1][x+1];
          const Real distPz = SDFLAB[z+1+1][y+1][x+1];
          const Real distMz = SDFLAB[z+1-1][y+1][x+1];
          // gradU
          const Real gradUX = inv2h*(distPx - distMx);
          const Real gradUY = inv2h*(distPy - distMy);
          const Real gradUZ = inv2h*(distPz - distMz);
          const Real gradUSq = gradUX*gradUX+gradUY*gradUY+gradUZ*gradUZ + EPS;

          const Real gradHX = (x == 0) ? 2.0*(-0.5*CHI[z][y][x+2]+2.0*CHI[z][y][x+1]-1.5*CHI[z][y][x]) : ( (x==Nx-1) ? 2.0*(1.5*CHI[z][y][x]-2.0*CHI[z][y][x-1]+0.5*CHI[z][y][x-2]) : (CHI[z][y][x+1]-CHI[z][y][x-1]));
          const Real gradHY = (y == 0) ? 2.0*(-0.5*CHI[z][y+2][x]+2.0*CHI[z][y+1][x]-1.5*CHI[z][y][x]) : ( (y==Ny-1) ? 2.0*(1.5*CHI[z][y][x]-2.0*CHI[z][y-1][x]+0.5*CHI[z][y-2][x]) : (CHI[z][y+1][x]-CHI[z][y-1][x]));
          const Real gradHZ = (z == 0) ? 2.0*(-0.5*CHI[z+2][y][x]+2.0*CHI[z+1][y][x]-1.5*CHI[z][y][x]) : ( (z==Nz-1) ? 2.0*(1.5*CHI[z][y][x]-2.0*CHI[z-1][y][x]+0.5*CHI[z-2][y][x]) : (CHI[z+1][y][x]-CHI[z-1][y][x]));

          if (gradHX*gradHX + gradHY*gradHY + gradHZ*gradHZ < 1e-12) continue;

          const Real numD = gradHX*gradUX + gradHY*gradUY + gradHZ*gradUZ;
          const Real Delta = fac1 * numD/gradUSq; //h^3 * Delta
          if (Delta>EPS) o->write(x,y,z,Delta,gradUX,gradUY,gradUZ);
      }
      o->allocate_surface();
    }
  }
};

}  // anonymous namespace


/// Compute chi-based center of mass for each obstacle.
static void kernelComputeGridCoM(SimulationData &sim)
{
  // TODO: Refactor to use only one omp parallel and only one MPI_Allreduce.
  for (const auto &obstacle : sim.obstacle_vector->getObstacleVector()) {
    Real com[4] = {0.0, 0.0, 0.0, 0.0};
    const auto& obstblocks = obstacle->getObstacleBlocks();
    #pragma omp parallel for schedule(static,1) reduction(+ : com[:4])
    for (size_t i=0; i<obstblocks.size(); i++) {
      if(obstblocks[i] == nullptr) continue;
      com[0] += obstblocks[i]->mass;
      com[1] += obstblocks[i]->CoM_x;
      com[2] += obstblocks[i]->CoM_y;
      com[3] += obstblocks[i]->CoM_z;
    }
    MPI_Allreduce(MPI_IN_PLACE, com, 4,MPI_Real,MPI_SUM, sim.comm);

    assert(com[0]>std::numeric_limits<Real>::epsilon());
    obstacle->centerOfMass[0] = com[1]/com[0];
    obstacle->centerOfMass[1] = com[2]/com[0];
    obstacle->centerOfMass[2] = com[3]/com[0];
  }
}

static void _kernelIntegrateUdefMomenta(SimulationData& sim, const BlockInfo& info)
{
  const int Nx = VectorBlock::sizeX;
  const int Ny = VectorBlock::sizeY;
  const int Nz = VectorBlock::sizeZ;
  for (const auto &obstacle : sim.obstacle_vector->getObstacleVector()) {
    const auto& obstblocks = obstacle->getObstacleBlocks();
    ObstacleBlock*const o = obstblocks[info.blockID];
    if (o == nullptr) continue;

    const std::array<Real,3> CM = obstacle->getCenterOfMass();
    //We use last momentum computed by this method to stabilize the computation
    //of the ang vel. This is because J is really tiny.
    const std::array<Real,3> oldCorrVel = {{
      obstacle->transVel_correction[0],
      obstacle->transVel_correction[1],
      obstacle->transVel_correction[2]
    }};

    const CHIMAT & __restrict__ CHI = o->chi;
    const UDEFMAT & __restrict__ UDEF = o->udef;
    Real &VV = o->V;
    Real &FX = o->FX, &FY = o->FY, &FZ = o->FZ;
    Real &TX = o->TX, &TY = o->TY, &TZ = o->TZ;
    Real &J0 = o->J0, &J1 = o->J1, &J2 = o->J2;
    Real &J3 = o->J3, &J4 = o->J4, &J5 = o->J5;
    VV = 0; FX = 0; FY = 0; FZ = 0; TX = 0; TY = 0; TZ = 0;
    J0 = 0; J1 = 0; J2 = 0; J3 = 0; J4 = 0; J5 = 0;

    for(int z=0; z<Nz; ++z)
    for(int y=0; y<Ny; ++y)
    for(int x=0; x<Nx; ++x)
    {
      if (CHI[z][y][x] <= 0) continue;
      Real p[3]; info.pos(p,x,y,z);
      const Real dv = info.h*info.h*info.h, X = CHI[z][y][x];
      p[0] -= CM[0];
      p[1] -= CM[1];
      p[2] -= CM[2];
      const Real dUs = UDEF[z][y][x][0] - oldCorrVel[0];
      const Real dVs = UDEF[z][y][x][1] - oldCorrVel[1];
      const Real dWs = UDEF[z][y][x][2] - oldCorrVel[2];
      VV += X * dv;
      FX += X * UDEF[z][y][x][0] * dv;
      FY += X * UDEF[z][y][x][1] * dv;
      FZ += X * UDEF[z][y][x][2] * dv;
      TX += X * ( p[1]*dWs - p[2]*dVs ) * dv;
      TY += X * ( p[2]*dUs - p[0]*dWs ) * dv;
      TZ += X * ( p[0]*dVs - p[1]*dUs ) * dv;
      J0 += X * ( p[1]*p[1]+p[2]*p[2] ) * dv; J3 -= X * p[0]*p[1] * dv;
      J1 += X * ( p[0]*p[0]+p[2]*p[2] ) * dv; J4 -= X * p[0]*p[2] * dv;
      J2 += X * ( p[0]*p[0]+p[1]*p[1] ) * dv; J5 -= X * p[1]*p[2] * dv;
    }
  }
}

/// Integrate momenta over the grid.
static void kernelIntegrateUdefMomenta(SimulationData& sim)
{
  const std::vector<cubism::BlockInfo>& chiInfo = sim.chiInfo();
  #pragma omp parallel for schedule(dynamic, 1)
  for (size_t i = 0; i < chiInfo.size(); ++i)
    _kernelIntegrateUdefMomenta(sim, chiInfo[i]);
}

/// Reduce momenta across blocks and MPI.
static void kernelAccumulateUdefMomenta(SimulationData& sim, bool justDebug = false)
{
  // TODO: Refactor to use only one omp parallel and one MPI_Allreduce.
  for (const auto &obst : sim.obstacle_vector->getObstacleVector()) {
    Real M[13] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    const auto& oBlock = obst->getObstacleBlocks();
    #pragma omp parallel for schedule(static,1) reduction(+ : M[:13])
    for (size_t i=0; i<oBlock.size(); i++) {
      if(oBlock[i] == nullptr) continue;
      M[ 0] += oBlock[i]->V ;
      M[ 1] += oBlock[i]->FX; M[ 2] += oBlock[i]->FY; M[ 3] += oBlock[i]->FZ;
      M[ 4] += oBlock[i]->TX; M[ 5] += oBlock[i]->TY; M[ 6] += oBlock[i]->TZ;
      M[ 7] += oBlock[i]->J0; M[ 8] += oBlock[i]->J1; M[ 9] += oBlock[i]->J2;
      M[10] += oBlock[i]->J3; M[11] += oBlock[i]->J4; M[12] += oBlock[i]->J5;
    }
    const auto comm = sim.comm;
    MPI_Allreduce(MPI_IN_PLACE, M, 13, MPI_Real, MPI_SUM, comm);
    assert(M[0] > EPS);

    const GenV AM = {{ M[ 4], M[ 5], M[ 6] }};
    const SymM J =  {{ M[ 7], M[ 8], M[ 9], M[10], M[11], M[12] }};
    const SymM invJ = invertSym(J);

    if(justDebug) {
      assert(std::fabs(M[ 1])<100*EPS);
      assert(std::fabs(M[ 2])<100*EPS);
      assert(std::fabs(M[ 3])<100*EPS);
      assert(std::fabs(AM[0])<100*EPS);
      assert(std::fabs(AM[1])<100*EPS);
      assert(std::fabs(AM[2])<100*EPS);
    } else {
      //solve avel = invJ \dot angMomentum
      obst->mass                   = M[ 0];
      obst->transVel_correction[0] = M[ 1] / M[0];
      obst->transVel_correction[1] = M[ 2] / M[0];
      obst->transVel_correction[2] = M[ 3] / M[0];
      obst->J[0] = M[ 7]; obst->J[1] = M[ 8]; obst->J[2] = M[ 9];
      obst->J[3] = M[10]; obst->J[4] = M[11]; obst->J[5] = M[12];
      obst->angVel_correction[0] = invJ[0]*AM[0] +invJ[3]*AM[1] +invJ[4]*AM[2];
      obst->angVel_correction[1] = invJ[3]*AM[0] +invJ[1]*AM[1] +invJ[5]*AM[2];
      obst->angVel_correction[2] = invJ[4]*AM[0] +invJ[5]*AM[1] +invJ[2]*AM[2];
    }
  }
}

/// Remove momenta from udef.
static void kernelRemoveUdefMomenta(SimulationData& sim, bool justDebug = false)
{
  const int Nx = VectorBlock::sizeX;
  const int Ny = VectorBlock::sizeY;
  const int Nz = VectorBlock::sizeZ;
  const std::vector<BlockInfo>& chiInfo = sim.chiInfo();

  // TODO: Refactor to use only one omp parallel.
  for (const auto &obstacle : sim.obstacle_vector->getObstacleVector()) {
    const std::array<Real, 3> angVel_correction = obstacle->angVel_correction;
    const std::array<Real, 3> transVel_correction = obstacle->transVel_correction;

    const std::array<Real,3> CM = obstacle->getCenterOfMass();
    const auto & obstacleBlocks = obstacle->getObstacleBlocks();

    #pragma omp parallel for schedule(dynamic, 1)
    for(size_t i=0; i < chiInfo.size(); i++)
    {
      const BlockInfo& info = chiInfo[i];
      const auto pos = obstacleBlocks[info.blockID];
      if(pos == nullptr) continue;
      UDEFMAT & __restrict__ UDEF = pos->udef;
      for(int z=0; z<Nz; ++z)
      for(int y=0; y<Ny; ++y)
      for(int x=0; x<Nx; ++x)
      {
        Real p[3]; info.pos(p,x,y,z);
        p[0] -= CM[0]; p[1] -= CM[1]; p[2] -= CM[2];
        const Real rotVel_correction[3] = {
          angVel_correction[1]*p[2] - angVel_correction[2]*p[1],
          angVel_correction[2]*p[0] - angVel_correction[0]*p[2],
          angVel_correction[0]*p[1] - angVel_correction[1]*p[0]
        };
        UDEF[z][y][x][0] -= transVel_correction[0] + rotVel_correction[0];
        UDEF[z][y][x][1] -= transVel_correction[1] + rotVel_correction[1];
        UDEF[z][y][x][2] -= transVel_correction[2] + rotVel_correction[2];
      }
    }
  }
}

void CreateObstacles::operator()(const Real dt)
{
  if(sim.obstacle_vector->nObstacles() == 0) return;
  if(sim.MeshChanged == false && sim.StaticObstacles) return;
  sim.MeshChanged = false;

  std::vector<BlockInfo>& chiInfo = sim.chiInfo();
  #pragma omp parallel for schedule(static)
  for (size_t i = 0; i < chiInfo.size(); ++i)
  {
    ScalarBlock& CHI = *(ScalarBlock*)chiInfo[i].ptrBlock;
    CHI.clear();
  }

  // Obstacles' advection must be done after we perform penalization:
  sim.uinf = sim.obstacle_vector->updateUinf();
  sim.obstacle_vector->update();

  { // put signed distance function on the grid
    sim.obstacle_vector->create();
  }

  {
    #pragma omp parallel
    {
      auto vecOB = sim.obstacle_vector->getAllObstacleBlocks();
      const KernelCharacteristicFunction K(vecOB);
      #pragma omp for
      for (size_t i = 0; i < chiInfo.size(); ++i)
      {
        ScalarBlock& CHI = *(ScalarBlock*)chiInfo[i].ptrBlock;
        K.operate(chiInfo[i],CHI);
      }
    }
  }

  // compute actual CoM given the CHI on the grid
  kernelComputeGridCoM(sim);
  kernelIntegrateUdefMomenta(sim);
  kernelAccumulateUdefMomenta(sim);
  kernelRemoveUdefMomenta(sim);
  sim.obstacle_vector->finalize(); // whatever else the obstacle needs
}

CubismUP_3D_NAMESPACE_END
