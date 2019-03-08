//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#include "operators/ObstaclesUpdate.h"
#include "obstacles/ObstacleVector.h"
#include "utils/MatArrayMath.h"

CubismUP_3D_NAMESPACE_BEGIN
using namespace cubism;

using CHIMAT = Real[CUP_BLOCK_SIZE][CUP_BLOCK_SIZE][CUP_BLOCK_SIZE];
using UDEFMAT = Real[CUP_BLOCK_SIZE][CUP_BLOCK_SIZE][CUP_BLOCK_SIZE][3];
static constexpr Real EPS = std::numeric_limits<Real>::epsilon();
static constexpr Real DBLEPS = std::numeric_limits<double>::epsilon();

struct KernelIntegrateFluidMomenta : public ObstacleVisitor
{
  const double lambda, dt;
  ObstacleVector * const obstacle_vector;
  const cubism::BlockInfo * info_ptr = nullptr;
  inline double dvol(const cubism::BlockInfo&info, const int x, const int y, const int z) const {
    double h[3]; info.spacing(h, x, y, z);
    return h[0] * h[1] * h[2];
  }

  KernelIntegrateFluidMomenta(double _dt, double _lambda, ObstacleVector* ov)
    : lambda(_lambda), dt(_dt), obstacle_vector(ov) {}

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

  void visit(Obstacle* const op)
  {
    const BlockInfo& info = * info_ptr;
    assert(info_ptr not_eq nullptr);
    const std::vector<ObstacleBlock*>& obstblocks = op->getObstacleBlocks();
    ObstacleBlock*const o = obstblocks[info.blockID];
    if (o == nullptr) return;

    const std::array<double,3> CM = op->getCenterOfMass();
    const FluidBlock &b = *(FluidBlock *)info.ptrBlock;
    const CHIMAT & __restrict__ CHI = o->chi;
    double &VV = o->V;
    double &FX = o->FX, &FY = o->FY, &FZ = o->FZ;
    double &TX = o->TX, &TY = o->TY, &TZ = o->TZ;
    VV = 0; FX = 0; FY = 0; FZ = 0; TX = 0; TY = 0; TZ = 0;
    double &J0 = o->J0, &J1 = o->J1, &J2 = o->J2;
    double &J3 = o->J3, &J4 = o->J4, &J5 = o->J5;
    J0 = 0; J1 = 0; J2 = 0; J3 = 0; J4 = 0; J5 = 0;

    for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
    for(int iy=0; iy<FluidBlock::sizeY; ++iy)
    for(int ix=0; ix<FluidBlock::sizeX; ++ix)
    {
      if (CHI[iz][iy][ix] <= 0) continue;
      double p[3]; info.pos(p, ix, iy, iz);
      const double dv = dvol(info, ix, iy, iz), X = CHI[iz][iy][ix];
      p[0] -= CM[0]; p[1] -= CM[1]; p[2] -= CM[2];

      VV += X * dv;
      J0 += X * dv * ( p[1]*p[1] + p[2]*p[2] );
      J1 += X * dv * ( p[0]*p[0] + p[2]*p[2] );
      J2 += X * dv * ( p[0]*p[0] + p[1]*p[1] );
      J3 -= X * dv * p[0]*p[1];
      J4 -= X * dv * p[0]*p[2];
      J5 -= X * dv * p[1]*p[2];

      FX += X * dv * b(ix,iy,iz).u;
      FY += X * dv * b(ix,iy,iz).v;
      FZ += X * dv * b(ix,iy,iz).w;
      TX += X * dv * ( p[1] * b(ix,iy,iz).w - p[2] * b(ix,iy,iz).v );
      TY += X * dv * ( p[2] * b(ix,iy,iz).u - p[0] * b(ix,iy,iz).w );
      TZ += X * dv * ( p[0] * b(ix,iy,iz).v - p[1] * b(ix,iy,iz).u );
    }
  }
};

struct KernelFinalizeObstacleVel : public ObstacleVisitor
{
  const double dt, lambda;
  FluidGridMPI * const grid;

  std::array<   int, 3> nSum = {{0, 0, 0}};
  std::array<double, 3> uSum = {{0, 0, 0}};

  KernelFinalizeObstacleVel(double _dt, double _lambda, FluidGridMPI*g) :
    dt(_dt), lambda(_lambda), grid(g) { }

  void visit(Obstacle* const obst)
  {
    double M[13] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
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
    assert(std::fabs(obst->mass - M[ 0]) < 10*DBLEPS);
    assert(std::fabs(obst->J[0] - M[ 7]) < 10*DBLEPS);
    assert(std::fabs(obst->J[1] - M[ 8]) < 10*DBLEPS);
    assert(std::fabs(obst->J[2] - M[ 9]) < 10*DBLEPS);
    assert(std::fabs(obst->J[3] - M[10]) < 10*DBLEPS);
    assert(std::fabs(obst->J[4] - M[11]) < 10*DBLEPS);
    assert(std::fabs(obst->J[5] - M[12]) < 10*DBLEPS);
    assert(M[0] > DBLEPS);
    const double mass = M[0];
    const int bRotate[3] = { obst->bBlockRotation[0] ? 0 : 1,
                             obst->bBlockRotation[1] ? 0 : 1,
                             obst->bBlockRotation[2] ? 0 : 1 };
    // Obstacles may prevent rotation along certain coordinates, but we want to
    // still record the imposed angular acceleration for dbg/post.
    // Cross terms in the J matrix describe how angular momentum is transferred
    // across directions. If direction is locked it absorbs energy (ie. applies
    // constraint force) which prevents momentum transfer from that direction.
    // We accomplish this by zeroing some terms in invJ. E.g if obstacle can
    // only rotate around z, \dot{omega}_z = \tau_z / J_zz
    // If obstacle does not rotate in any direction we still want to see accels
    // (for debug) and therefore diagonal is never zeroed.
    const GenM fullJ = {{
      M[ 7]             , M[10] * bRotate[1], M[11] * bRotate[2],
      M[10] * bRotate[0], M[ 8]             , M[12] * bRotate[2],
      M[11] * bRotate[0], M[12] * bRotate[1], M[ 9]
    }};
    const GenM invJ = invertGen(fullJ);
    if(bRotate[0]==0 && bRotate[1]==0) // J_zz * invJzz == 1
      assert(std::fabs(invJ[8]*M[9]-1) < 10*DBLEPS);

    // these are the angular and linear momenta of the fluid inside the obstacle
    obst->angVel_fluid = multGenVec(invJ, GenV{{ M[4], M[5], M[6] }});
    obst->transVel_fluid = {M[1] / mass, M[2] / mass, M[3] / mass};

    #ifndef OLD_INTEGRATE_MOM
      const Real penalFac = lambda*dt>32? 1 : 1 - std::exp(- lambda * dt);
      // compute linear and angular impulses onto obstacle:
      const Real impLX = penalFac*(obst->transVel_fluid[0] - obst->transVel[0]);
      const Real impLY = penalFac*(obst->transVel_fluid[1] - obst->transVel[1]);
      const Real impLZ = penalFac*(obst->transVel_fluid[2] - obst->transVel[2]);
      const Real impAX = penalFac*(obst->angVel_fluid[0]   - obst->angVel[0]  );
      const Real impAY = penalFac*(obst->angVel_fluid[1]   - obst->angVel[1]  );
      const Real impAZ = penalFac*(obst->angVel_fluid[2]   - obst->angVel[2]  );
      obst->force[0]  = impLX * mass / dt;
      obst->force[1]  = impLY * mass / dt;
      obst->force[2]  = impLZ * mass / dt;
      obst->torque[0] = (fullJ[0]*impAX + fullJ[1]*impAY + fullJ[2]*impAZ)/dt;
      obst->torque[1] = (fullJ[3]*impAX + fullJ[4]*impAY + fullJ[5]*impAZ)/dt;
      obst->torque[2] = (fullJ[6]*impAX + fullJ[7]*impAY + fullJ[8]*impAZ)/dt;
      obst->transVel_computed[0] = obst->transVel[0] + impLX;
      obst->transVel_computed[1] = obst->transVel[1] + impLY;
      obst->transVel_computed[2] = obst->transVel[2] + impLZ;
      obst->angVel_computed[0]   = obst->angVel[0]   + impAX;
      obst->angVel_computed[1]   = obst->angVel[1]   + impAY;
      obst->angVel_computed[2]   = obst->angVel[2]   + impAZ;
    #else // momenta in the flow pre penal become momenta of the obstacle
      obst->force[0]  = M[1]; //these are all wrong now, must compute post penal
      obst->force[1]  = M[2]; //these are all wrong now, must compute post penal
      obst->force[2]  = M[3]; //these are all wrong now, must compute post penal
      obst->torque[0] = M[4]; //these are all wrong now, must compute post penal
      obst->torque[1] = M[5]; //these are all wrong now, must compute post penal
      obst->torque[2] = M[6]; //these are all wrong now, must compute post penal
      obst->transVel_computed[0] = obst->transVel_fluid[0];
      obst->transVel_computed[1] = obst->transVel_fluid[1];
      obst->transVel_computed[2] = obst->transVel_fluid[2];
      obst->angVel_computed[0] = obst->angVel_fluid[0];
      obst->angVel_computed[1] = obst->angVel_fluid[1];
      obst->angVel_computed[2] = obst->angVel_fluid[2];
    #endif

    obst->computeVelocities();

    const auto &bFixFrameOfRef = obst->bFixFrameOfRef;
    if (bFixFrameOfRef[0]) { nSum[0]++; uSum[0] -= obst->transVel[0]; }
    if (bFixFrameOfRef[1]) { nSum[1]++; uSum[1] -= obst->transVel[1]; }
    if (bFixFrameOfRef[2]) { nSum[2]++; uSum[2] -= obst->transVel[2]; }
  }
};

void UpdateObstacles::operator()(const double dt)
{
  if(sim.obstacle_vector->nObstacles() == 0) return;

  sim.startProfiler("Obst Int Vel");
  { // integrate momenta by looping over grid
    #pragma omp parallel
    { // each thread needs to call its own non-const operator() function
      auto K = KernelIntegrateFluidMomenta(dt, sim.lambda, sim.obstacle_vector);
      #pragma omp for schedule(dynamic, 1)
      for (size_t i = 0; i < vInfo.size(); ++i) K(vInfo[i]);
    }
  }
  sim.stopProfiler();

  sim.startProfiler("Obst Upd Vel");
  auto* K = new KernelFinalizeObstacleVel(dt, sim.lambda, sim.grid);
  ObstacleVisitor* kernel = static_cast<ObstacleVisitor*>(K);
  assert(kernel not_eq nullptr);
  sim.obstacle_vector->Accept(kernel); // accept you son of a french cow
  if( K->nSum[0] > 0 ) sim.uinf[0] = K->uSum[0] / K->nSum[0];
  if( K->nSum[1] > 0 ) sim.uinf[1] = K->uSum[1] / K->nSum[1];
  if( K->nSum[2] > 0 ) sim.uinf[2] = K->uSum[2] / K->nSum[2];
  //if(rank == 0) if(nSum[0] || nSum[1] || nSum[2])
  //  printf("New Uinf %g %g %g (from %d %d %d)\n",
  //  uInf[0],uInf[1],uInf[2],nSum[0],nSum[1],nSum[2]);
  delete K;
  // Obstacles' advection must be done after we compute new uinf :
  sim.obstacle_vector->update();
  sim.stopProfiler();

  check("UpdateObstacles");
}

CubismUP_3D_NAMESPACE_END

#if 0
const double UFROT[3] = {
    fluidAngVel[1]*p[2] - fluidAngVel[2]*p[1],
    fluidAngVel[2]*p[0] - fluidAngVel[0]*p[2],
    fluidAngVel[0]*p[1] - fluidAngVel[1]*p[0]
};
const double USROT[3] = {
    omega[1]*p[2] - omega[2]*p[1],
    omega[2]*p[0] - omega[0]*p[2],
    omega[0]*p[1] - omega[1]*p[0]
};
const Real velFluctX = b(ix,iy,iz).u - fluidLinVel[0] - UFROT[0];
const Real velFluctY = b(ix,iy,iz).v - fluidLinVel[1] - UFROT[1];
const Real velFluctZ = b(ix,iy,iz).w - fluidLinVel[2] - UFROT[2];
const Real deltaLinX = (            UT[0]-fluidLinVel[0]) * penalFX;
const Real deltaLinY = (            UT[1]-fluidLinVel[1]) * penalFY;
const Real deltaLinZ = (            UT[2]-fluidLinVel[2]) * penalFZ;
const Real deltaAngX = (         USROT[0]-      UFROT[0]) * penalTX;
const Real deltaAngY = (         USROT[1]-      UFROT[1]) * penalTY;
const Real deltaAngZ = (         USROT[2]-      UFROT[2]) * penalTZ;
const Real deltaDefX = (UD[iz][iy][ix][0]-    velFluctX ) * penalDef;
const Real deltaDefY = (UD[iz][iy][ix][1]-    velFluctY ) * penalDef;
const Real deltaDefZ = (UD[iz][iy][ix][2]-    velFluctZ ) * penalDef;
b(ix,iy,iz).u += CHI[iz][iy][ix]*(deltaLinX + deltaAngX + deltaDefX);
b(ix,iy,iz).v += CHI[iz][iy][ix]*(deltaLinY + deltaAngY + deltaDefY);
b(ix,iy,iz).w += CHI[iz][iy][ix]*(deltaLinZ + deltaAngZ + deltaDefZ);
#endif
