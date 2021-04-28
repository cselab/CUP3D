//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#include "ObstaclesCreate.h"
#include "../obstacles/ObstacleVector.h"
#include "../utils/MatArrayMath.h"

CubismUP_3D_NAMESPACE_BEGIN
using namespace cubism;

namespace {

using CHIMAT =  Real[CUP_BLOCK_SIZEZ][CUP_BLOCK_SIZEY][CUP_BLOCK_SIZEX];
using UDEFMAT = Real[CUP_BLOCK_SIZEZ][CUP_BLOCK_SIZEY][CUP_BLOCK_SIZEX][3];
static constexpr Real EPS = std::numeric_limits<Real>::epsilon();

class KernelCharacteristicFunction
{
  using v_v_ob = std::vector<std::vector<ObstacleBlock*>*>;
  const v_v_ob & vec_obstacleBlocks;

  public:
  KernelCharacteristicFunction(const v_v_ob& v) : vec_obstacleBlocks(v) {}

  Real J(const Real x) const
  {
    return (x > 0.0) ? 0.5*x*x : 0.0;
  }
  Real I(const Real x) const
  {
    return std::max(x,0.0);
  }

  template <typename BlockType>
  void operate(const BlockInfo& info, BlockType& b) const
  {
    const Real h = info.h_gridpoint, inv2h = .5/h, fac1 = .5*h*h, vol = h*h*h;
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
      for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
      for(int iy=0; iy<FluidBlock::sizeY; ++iy)
      for(int ix=0; ix<FluidBlock::sizeX; ++ix)
      {
        #if 1
        // here I read fist from SDF to deal with obstacles sharing block
        if (SDFLAB[iz+1][iy+1][ix+1] > +gp*h || SDFLAB[iz+1][iy+1][ix+1] < -gp*h)
        {
          CHI[iz][iy][ix] = SDFLAB[iz+1][iy+1][ix+1] > 0 ? 1 : 0;
        }
        else
        {
          const Real distPx = SDFLAB[iz+1][iy+1][ix+1+1];
          const Real distMx = SDFLAB[iz+1][iy+1][ix+1-1];
          const Real distPy = SDFLAB[iz+1][iy+1+1][ix+1];
          const Real distMy = SDFLAB[iz+1][iy+1-1][ix+1];
          const Real distPz = SDFLAB[iz+1+1][iy+1][ix+1];
          const Real distMz = SDFLAB[iz+1-1][iy+1][ix+1];
          // gradU
          const Real gradUX = inv2h*(distPx - distMx);
          const Real gradUY = inv2h*(distPy - distMy);
          const Real gradUZ = inv2h*(distPz - distMz);
          const Real gradUSq = gradUX*gradUX+gradUY*gradUY+gradUZ*gradUZ + EPS;
          const Real IplusX = std::max(0.0,distPx);
          const Real IminuX = std::max(0.0,distMx);
          const Real IplusY = std::max(0.0,distPy);
          const Real IminuY = std::max(0.0,distMy);
          const Real IplusZ = std::max(0.0,distPz);
          const Real IminuZ = std::max(0.0,distMz);
          // gradI: first primitive of H(x): I(x) = int_0^x H(y) dy
          const Real gradIX = inv2h*(IplusX - IminuX);
          const Real gradIY = inv2h*(IplusY - IminuY);
          const Real gradIZ = inv2h*(IplusZ - IminuZ);
          const Real numH = gradIX*gradUX + gradIY*gradUY + gradIZ*gradUZ;
          CHI[iz][iy][ix] = numH/gradUSq;
        }
        #else
          CHI[iz][iy][ix] = SDFLAB[iz+1][iy+1][ix+1] > 0 ? 1 : 0;
        #endif
        Real p[3]; info.pos(p, ix,iy,iz);
        b(ix,iy,iz).chi = std::max(CHI[iz][iy][ix], b(ix,iy,iz).chi);
        o->CoM_x += CHI[iz][iy][ix] * vol * p[0];
        o->CoM_y += CHI[iz][iy][ix] * vol * p[1];
        o->CoM_z += CHI[iz][iy][ix] * vol * p[2];
        o->mass  += CHI[iz][iy][ix] * vol;
      }

      for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
      for(int iy=0; iy<FluidBlock::sizeY; ++iy)
      for(int ix=0; ix<FluidBlock::sizeX; ++ix)
      {
          const Real distPx = SDFLAB[iz+1][iy+1][ix+1+1];
          const Real distMx = SDFLAB[iz+1][iy+1][ix+1-1];
          const Real distPy = SDFLAB[iz+1][iy+1+1][ix+1];
          const Real distMy = SDFLAB[iz+1][iy+1-1][ix+1];
          const Real distPz = SDFLAB[iz+1+1][iy+1][ix+1];
          const Real distMz = SDFLAB[iz+1-1][iy+1][ix+1];
          // gradU
          const Real gradUX = inv2h*(distPx - distMx);
          const Real gradUY = inv2h*(distPy - distMy);
          const Real gradUZ = inv2h*(distPz - distMz);
          const Real gradUSq = gradUX*gradUX+gradUY*gradUY+gradUZ*gradUZ + EPS;

          const Real gradHX = (ix == 0) ? 2.0*(-0.5*CHI[iz][iy][ix+2]+2.0*CHI[iz][iy][ix+1]-1.5*CHI[iz][iy][ix]) : ( (ix==FluidBlock::sizeX-1) ? 2.0*(1.5*CHI[iz][iy][ix]-2.0*CHI[iz][iy][ix-1]+0.5*CHI[iz][iy][ix-2]) : (CHI[iz][iy][ix+1]-CHI[iz][iy][ix-1]));
          const Real gradHY = (iy == 0) ? 2.0*(-0.5*CHI[iz][iy+2][ix]+2.0*CHI[iz][iy+1][ix]-1.5*CHI[iz][iy][ix]) : ( (iy==FluidBlock::sizeY-1) ? 2.0*(1.5*CHI[iz][iy][ix]-2.0*CHI[iz][iy-1][ix]+0.5*CHI[iz][iy-2][ix]) : (CHI[iz][iy+1][ix]-CHI[iz][iy-1][ix]));
          const Real gradHZ = (iz == 0) ? 2.0*(-0.5*CHI[iz+2][iy][ix]+2.0*CHI[iz+1][iy][ix]-1.5*CHI[iz][iy][ix]) : ( (iz==FluidBlock::sizeZ-1) ? 2.0*(1.5*CHI[iz][iy][ix]-2.0*CHI[iz-1][iy][ix]+0.5*CHI[iz-2][iy][ix]) : (CHI[iz+1][iy][ix]-CHI[iz-1][iy][ix]));

          const Real numD = gradHX*gradUX + gradHY*gradUY + gradHZ*gradUZ;
          const Real Delta = fac1 * numD/gradUSq; //h^3 * Delta
          if (Delta>EPS) o->write(ix, iy, iz, Delta, gradUX, gradUY, gradUZ);
      }
      o->allocate_surface();
    }
  }
};

struct KernelComputeGridCoM : public ObstacleVisitor
{
  FluidGridMPI * const grid;
  const std::vector<cubism::BlockInfo>& vInfo = grid->getBlocksInfo();

  KernelComputeGridCoM(FluidGridMPI*g) : grid(g) { }

  void visit(Obstacle* const obstacle)
  {
    double com[4] = {0.0, 0.0, 0.0, 0.0};
    const auto& obstblocks = obstacle->getObstacleBlocks();
    #pragma omp parallel for schedule(static,1) reduction(+ : com[:4])
    for (size_t i=0; i<obstblocks.size(); i++) {
      if(obstblocks[i] == nullptr) continue;
      com[0] += obstblocks[i]->mass;
      com[1] += obstblocks[i]->CoM_x;
      com[2] += obstblocks[i]->CoM_y;
      com[3] += obstblocks[i]->CoM_z;
    }
    MPI_Allreduce(MPI_IN_PLACE, com, 4,MPI_DOUBLE,MPI_SUM, grid->getCartComm());

    assert(com[0]>std::numeric_limits<Real>::epsilon());
    obstacle->centerOfMass[0] = com[1]/com[0];
    obstacle->centerOfMass[1] = com[2]/com[0];
    obstacle->centerOfMass[2] = com[3]/com[0];
  }
};

struct KernelIntegrateUdefMomenta : public ObstacleVisitor
{
  ObstacleVector * const obstacle_vector;
  const cubism::BlockInfo * info_ptr = nullptr;
  inline double dvol(const cubism::BlockInfo&info, const int x, const int y, const int z) const {
    double h[3]; info.spacing(h, x, y, z);
    return h[0] * h[1] * h[2];
  }

  KernelIntegrateUdefMomenta(ObstacleVector* ov) : obstacle_vector(ov) {}

  void operator()(const cubism::BlockInfo& info)
  {
    // first store the lab and info, then do visitor
    assert(info_ptr == nullptr);
    info_ptr = & info;
    ObstacleVisitor* const base = static_cast<ObstacleVisitor*> (this);
    assert( base not_eq nullptr );
    obstacle_vector->Accept( base );
    info_ptr = nullptr;
  }

  void visit(Obstacle* const obstacle)
  {
    const BlockInfo& info = * info_ptr;
    assert(info_ptr not_eq nullptr);
    const auto& obstblocks = obstacle->getObstacleBlocks();
    ObstacleBlock*const o = obstblocks[info.blockID];
    if (o == nullptr) return;

    const std::array<double,3> CM = obstacle->getCenterOfMass();
    //We use last momentum computed by this method to stabilize the computation
    //of the ang vel. This is because J is really tiny.
    const std::array<double,3> oldCorrVel = {{
      obstacle->transVel_correction[0],
      obstacle->transVel_correction[1],
      obstacle->transVel_correction[2]
    }};

    const CHIMAT & __restrict__ CHI = o->chi;
    const UDEFMAT & __restrict__ UDEF = o->udef;
    double &VV = o->V;
    double &FX = o->FX, &FY = o->FY, &FZ = o->FZ;
    double &TX = o->TX, &TY = o->TY, &TZ = o->TZ;
    double &J0 = o->J0, &J1 = o->J1, &J2 = o->J2;
    double &J3 = o->J3, &J4 = o->J4, &J5 = o->J5;
    VV = 0; FX = 0; FY = 0; FZ = 0; TX = 0; TY = 0; TZ = 0;
    J0 = 0; J1 = 0; J2 = 0; J3 = 0; J4 = 0; J5 = 0;

    for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
    for(int iy=0; iy<FluidBlock::sizeY; ++iy)
    for(int ix=0; ix<FluidBlock::sizeX; ++ix)
    {
      if (CHI[iz][iy][ix] <= 0) continue;
      double p[3]; info.pos(p, ix, iy, iz);
      const double dv = dvol(info, ix, iy, iz), X = CHI[iz][iy][ix];
      p[0] -= CM[0];
      p[1] -= CM[1];
      p[2] -= CM[2];
      const double dUs = UDEF[iz][iy][ix][0] - oldCorrVel[0];
      const double dVs = UDEF[iz][iy][ix][1] - oldCorrVel[1];
      const double dWs = UDEF[iz][iy][ix][2] - oldCorrVel[2];
      VV += X * dv;
      FX += X * UDEF[iz][iy][ix][0] * dv;
      FY += X * UDEF[iz][iy][ix][1] * dv;
      FZ += X * UDEF[iz][iy][ix][2] * dv;
      TX += X * ( p[1]*dWs - p[2]*dVs ) * dv;
      TY += X * ( p[2]*dUs - p[0]*dWs ) * dv;
      TZ += X * ( p[0]*dVs - p[1]*dUs ) * dv;
      J0 += X * ( p[1]*p[1]+p[2]*p[2] ) * dv; J3 -= X * p[0]*p[1] * dv;
      J1 += X * ( p[0]*p[0]+p[2]*p[2] ) * dv; J4 -= X * p[0]*p[2] * dv;
      J2 += X * ( p[0]*p[0]+p[1]*p[1] ) * dv; J5 -= X * p[1]*p[2] * dv;
    }
  }
};

struct KernelAccumulateUdefMomenta : public ObstacleVisitor
{
  FluidGridMPI * const grid;
  const std::vector<cubism::BlockInfo>& vInfo = grid->getBlocksInfo();
  const bool justDebug = false;

  KernelAccumulateUdefMomenta(FluidGridMPI*g, bool dbg = false) :
    grid(g), justDebug(dbg) {}

  void visit(Obstacle* const obst)
  {
    double M[13] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
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
    const auto comm = grid->getCartComm();
    MPI_Allreduce(MPI_IN_PLACE, M, 13, MPI_DOUBLE, MPI_SUM, comm);
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
};

struct KernelRemoveUdefMomenta : public ObstacleVisitor
{
  FluidGridMPI * const grid;
  const std::vector<cubism::BlockInfo>& vInfo = grid->getBlocksInfo();

  KernelRemoveUdefMomenta(FluidGridMPI*g) : grid(g) { }

  void visit(Obstacle* const obstacle)
  {
    const std::array<double,3> angVel_correction = {{
      obstacle->angVel_correction[0],
      obstacle->angVel_correction[1],
      obstacle->angVel_correction[2]
    }};
    const std::array<double,3> transVel_correction = {{
      obstacle->transVel_correction[0],
      obstacle->transVel_correction[1],
      obstacle->transVel_correction[2]
    }};

    #ifdef CUP_VERBOSE
     if(sim.rank==0)
        printf("Obstacle %d moment corrections lin:[%f %f %f] ang:[%f %f %f]\n",
          obstacle->obstacleID, transVel_correction[0], transVel_correction[1],
          transVel_correction[2], angVel_correction[0], angVel_correction[1],
          angVel_correction[2]);
    #endif

    const std::array<double,3> CM = obstacle->getCenterOfMass();
    const auto & obstacleBlocks = obstacle->getObstacleBlocks();

    #pragma omp parallel for schedule(dynamic, 1)
    for(size_t i=0; i < vInfo.size(); i++)
    {
      const BlockInfo& info = vInfo[i];
      const auto pos = obstacleBlocks[info.blockID];
      if(pos == nullptr) continue;
      UDEFMAT & __restrict__ UDEF = pos->udef;
      for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
      for(int iy=0; iy<FluidBlock::sizeY; ++iy)
      for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
        double p[3]; info.pos(p, ix, iy, iz);
        p[0] -= CM[0]; p[1] -= CM[1]; p[2] -= CM[2];
        const double rotVel_correction[3] = {
          angVel_correction[1]*p[2] - angVel_correction[2]*p[1],
          angVel_correction[2]*p[0] - angVel_correction[0]*p[2],
          angVel_correction[0]*p[1] - angVel_correction[1]*p[0]
        };
        UDEF[iz][iy][ix][0] -= transVel_correction[0] + rotVel_correction[0];
        UDEF[iz][iy][ix][1] -= transVel_correction[1] + rotVel_correction[1];
        UDEF[iz][iy][ix][2] -= transVel_correction[2] + rotVel_correction[2];
      }
    }
  }
};

}

void CreateObstacles::operator()(const double dt)
{
  if(sim.obstacle_vector->nObstacles() == 0) return;

  sim.startProfiler("Obst Reset");
  std::vector<cubism::BlockInfo>& vInfo = sim.vInfo();
  #pragma omp parallel for schedule(static)
  for (size_t i = 0; i < vInfo.size(); ++i)
  {
    FluidBlock& b = *(FluidBlock*)vInfo[i].ptrBlock;
    for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
    for(int iy=0; iy<FluidBlock::sizeY; ++iy)
    for(int ix=0; ix<FluidBlock::sizeX; ++ix)
      b(ix,iy,iz).chi = 0; // will be accessed by max with pos def qtity
  }
  sim.stopProfiler();

  // Obstacles' advection must be done after we perform penalization:
  sim.uinf = sim.obstacle_vector->updateUinf();
  sim.obstacle_vector->update();

  sim.startProfiler("Obst SDF");
  { // put signed distance function on the grid
    sim.obstacle_vector->create();
  }
  sim.stopProfiler();

  sim.startProfiler("Obst CHI");
  {
    #pragma omp parallel
    {
      auto vecOB = sim.obstacle_vector->getAllObstacleBlocks();
      const KernelCharacteristicFunction K(vecOB);
      #pragma omp for
      for (size_t i = 0; i < vInfo.size(); ++i)
      {
        FluidBlock& b = *(FluidBlock*)vInfo[i].ptrBlock;
        K.operate<FluidBlock>(vInfo[i],b);
      }
    }
  }
  sim.stopProfiler();

  sim.startProfiler("Obst CoM");
  { // compute actual CoM given the CHI on the grid
    ObstacleVisitor* visitor = new KernelComputeGridCoM(grid);
    sim.obstacle_vector->Accept(visitor);
    delete visitor;
  }
  sim.stopProfiler();

  sim.startProfiler("Obst Int Mom");
  { // integrate momenta by looping over grid
    #pragma omp parallel
    { // each thread needs to call its own non-const operator() function
      KernelIntegrateUdefMomenta K(sim.obstacle_vector);
      #pragma omp for schedule(dynamic, 1)
      for (size_t i = 0; i < vInfo.size(); ++i) K(vInfo[i]);
    }
  }
  sim.stopProfiler();

  sim.startProfiler("Obst Rdx Mom");
  { // reduce momenta across blocks and MPI
    ObstacleVisitor* visitor = new KernelAccumulateUdefMomenta(grid);
    sim.obstacle_vector->Accept(visitor);
    delete visitor;
  }
  sim.stopProfiler();

  sim.startProfiler("Obst 0 Mom");
  { // remove momenta from udef
    ObstacleVisitor* visitor = new KernelRemoveUdefMomenta(grid);
    sim.obstacle_vector->Accept(visitor);
    delete visitor;
  }
  sim.stopProfiler();

  sim.startProfiler("Obst finalize");
  sim.obstacle_vector->finalize(); // whatever else the obstacle needs
  sim.stopProfiler();
  check("CreateObstacles");
}

CubismUP_3D_NAMESPACE_END
