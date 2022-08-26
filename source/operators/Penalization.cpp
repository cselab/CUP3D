//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//

#include "Penalization.h"
#include "../Obstacles/ObstacleVector.h"

CubismUP_3D_NAMESPACE_BEGIN
using namespace cubism;

namespace {

using CHIMAT =  Real[CUP_BLOCK_SIZEZ][CUP_BLOCK_SIZEY][CUP_BLOCK_SIZEX];
using UDEFMAT = Real[CUP_BLOCK_SIZEZ][CUP_BLOCK_SIZEY][CUP_BLOCK_SIZEX][3];

template<bool implicitPenalization>
struct KernelPenalization
{
  const Real dt, invdt = 1.0/dt, lambda;
  ObstacleVector * const obstacle_vector;

  KernelPenalization(Real _dt, Real _lambda, ObstacleVector* ov) :
    dt(_dt), lambda(_lambda), obstacle_vector(ov) {}

  void operator()(const cubism::BlockInfo& info,const BlockInfo& ChiInfo) const
  {
    for (const auto &obstacle : obstacle_vector->getObstacleVector())
      visit(info, ChiInfo, obstacle.get());
  }

  void visit(const BlockInfo& info, const BlockInfo& ChiInfo, Obstacle* const obstacle) const
  {
    const auto& obstblocks = obstacle->getObstacleBlocks();
    ObstacleBlock*const o = obstblocks[info.blockID];
    if (o == nullptr) return;

    const CHIMAT & __restrict__ CHI = o->chi;
    const UDEFMAT & __restrict__ UDEF = o->udef;
    VectorBlock& b = *(VectorBlock*)info.ptrBlock;
    ScalarBlock& bChi = *(ScalarBlock*)ChiInfo.ptrBlock;
    const std::array<Real,3> CM = obstacle->getCenterOfMass();
    const std::array<Real,3> vel = obstacle->getTranslationVelocity();
    const std::array<Real,3> omega = obstacle->getAngularVelocity();
    const Real dv = std::pow(info.h, 3);

    // Obstacle-specific lambda, useful for gradually adding an obstacle to the flow.
    // lambda = 1/dt hardcoded for expl time int, other options are wrong.
    const Real lambdaFac = implicitPenalization? lambda : invdt;

    Real &FX = o->FX, &FY = o->FY, &FZ = o->FZ;
    Real &TX = o->TX, &TY = o->TY, &TZ = o->TZ;
    FX = 0; FY = 0; FZ = 0; TX = 0; TY = 0; TZ = 0;

    for(int iz=0; iz<VectorBlock::sizeZ; ++iz)
    for(int iy=0; iy<VectorBlock::sizeY; ++iy)
    for(int ix=0; ix<VectorBlock::sizeX; ++ix)
    {
      // What if multiple obstacles share a block? Do not write udef onto
      // grid if CHI stored on the grid is greater than obst's CHI.
      if(bChi(ix,iy,iz).s > CHI[iz][iy][ix]) continue;
      if(CHI[iz][iy][ix] <= 0) continue; // no need to do anything
      Real p[3]; info.pos(p, ix, iy, iz);
      p[0] -= CM[0]; p[1] -= CM[1]; p[2] -= CM[2];

      const Real U_TOT[3] = {
          vel[0] + omega[1]*p[2] - omega[2]*p[1] + UDEF[iz][iy][ix][0],
          vel[1] + omega[2]*p[0] - omega[0]*p[2] + UDEF[iz][iy][ix][1],
          vel[2] + omega[0]*p[1] - omega[1]*p[0] + UDEF[iz][iy][ix][2]
      };
      const Real X = CHI[iz][iy][ix];
      const Real penalFac = implicitPenalization? X*lambdaFac/(1+X*lambdaFac*dt):X*lambdaFac;
      const Real FPX = penalFac * (U_TOT[0] - b(ix,iy,iz).u[0]);
      const Real FPY = penalFac * (U_TOT[1] - b(ix,iy,iz).u[1]);
      const Real FPZ = penalFac * (U_TOT[2] - b(ix,iy,iz).u[2]);
      // What if two obstacles overlap? Let's plus equal. We will need a
      // repulsion term of the velocity at some point in the code.
      b(ix,iy,iz).u[0] = b(ix,iy,iz).u[0] + dt * FPX;
      b(ix,iy,iz).u[1] = b(ix,iy,iz).u[1] + dt * FPY;
      b(ix,iy,iz).u[2] = b(ix,iy,iz).u[2] + dt * FPZ;

      FX += dv * FPX; FY += dv * FPY; FZ += dv * FPZ;
      TX += dv * ( p[1] * FPZ - p[2] * FPY );
      TY += dv * ( p[2] * FPX - p[0] * FPZ );
      TZ += dv * ( p[0] * FPY - p[1] * FPX );
    }
  }
};

static void kernelFinalizePenalizationForce(SimulationData& sim)
{
  // TODO: Refactor to use only one omp parallel and MPI_Allreduce.
  for (const auto &obst : sim.obstacle_vector->getObstacleVector()) {
    static constexpr int nQoI = 6;
    Real M[nQoI] = { 0 };
    const auto& oBlock = obst->getObstacleBlocks();
    #pragma omp parallel for schedule(static) reduction(+ : M[:nQoI])
    for (size_t i=0; i<oBlock.size(); ++i) {
      if(oBlock[i] == nullptr) continue;
      M[0] += oBlock[i]->FX; M[1] += oBlock[i]->FY; M[2] += oBlock[i]->FZ;
      M[3] += oBlock[i]->TX; M[4] += oBlock[i]->TY; M[5] += oBlock[i]->TZ;
    }
    const auto comm = sim.comm;
    MPI_Allreduce(MPI_IN_PLACE, M, nQoI, MPI_Real, MPI_SUM, comm);
    obst->force[0]  = M[0]; obst->force[1]  = M[1]; obst->force[2]  = M[2];
    obst->torque[0] = M[3]; obst->torque[1] = M[4]; obst->torque[2] = M[5];
  }
}

void ComputeJ(const Real * Rc, const Real * R, const Real * N, const Real * I, Real *J)
{
  //Computation of I^{-1} [(Rc - R) x N]

  //1. v = [(Rc - R) x N]
  const Real v0 = ( Rc[1] - R[1] )*N[2] - ( Rc[2] - R[2] )*N[1];
  const Real v1 = ( Rc[2] - R[2] )*N[0] - ( Rc[0] - R[0] )*N[2];
  const Real v2 = ( Rc[0] - R[0] )*N[1] - ( Rc[1] - R[1] )*N[0];

  //2. Invert I. Note that I[6] = {Ixx,Iyy,Izz,Ixy,Ixz,Iyz}
  const Real m00 = I[0];
  const Real m01 = I[3];
  const Real m02 = I[4];
  const Real m11 = I[1];
  const Real m12 = I[5];
  const Real m22 = I[2];

  //Invert the symmetric 3x3 matrix I. Inverse matrix a_ij is NOT multiplied by 1/det(I)
  const Real a00 = m22*m11 - m12*m12;
  const Real a01 = m02*m12 - m22*m01;
  const Real a02 = m01*m12 - m02*m11;
  const Real a11 = m22*m00 - m02*m02;
  const Real a12 = m01*m02 - m00*m12;
  const Real a22 = m00*m11 - m01*m01;

  //multiply final result with 1/det(I)
  const Real determinant =  1.0/((m00 * a00) + (m01 * a01) + (m02 * a02));

  J[0] = (a00*v0 + a01*v1 + a02*v2)*determinant;
  J[1] = (a01*v0 + a11*v1 + a12*v2)*determinant;
  J[2] = (a02*v0 + a12*v1 + a22*v2)*determinant;
}


void ElasticCollision(const Real m1, const Real m2, const Real*I1, const Real*I2, const Real *C1, const Real*C2,
                      const Real NX, const Real NY, const Real NZ, const Real CX, const Real  CY, const Real CZ,
                      const Real *vc1, const Real *vc2, Real *v1, Real *v2,Real *o1, Real *o2)
{
    const Real e = 1.0; // coefficient of restitution
    const Real N[3] ={NX,NY,NZ};
    const Real C[3] ={CX,CY,CZ};

    Real J1[3];
    Real J2[3]; 
    ComputeJ(C,C1,N,I1,J1);
    ComputeJ(C,C2,N,I2,J2);

    const Real nom = (e+1)*( (vc2[0]-vc1[0])*N[0] + (vc2[1]-vc1[1])*N[1] + (vc2[2]-vc1[2])*N[2] ); 

    const Real denom = (1.0/m1+1.0/m2) + 
                ( ( (J1[1]*(C[2]-C1[2])-J1[2]*(C[1]-C1[1])) + (J2[1]*(C[2]-C2[2])-J2[2]*(C[1]-C2[1]))) *N[0]+
                  ( (J1[2]*(C[0]-C1[0])-J1[0]*(C[2]-C1[2])) + (J2[2]*(C[0]-C2[0])-J2[0]*(C[2]-C2[2]))) *N[1]+
                  ( (J1[0]*(C[1]-C1[1])-J1[1]*(C[0]-C1[0])) + (J2[0]*(C[1]-C2[1])-J2[1]*(C[0]-C2[0]))) *N[2]);

    const Real impulse = nom/(denom+1e-21);

    v1[0] += (N[0]/m1)*impulse;
    v1[1] += (N[1]/m1)*impulse;
    v1[2] += (N[2]/m1)*impulse;

    v2[0] -= (N[0]/m2)*impulse;
    v2[1] -= (N[1]/m2)*impulse;
    v2[2] -= (N[2]/m2)*impulse;

    o1[0] += J1[0]*impulse;
    o1[1] += J1[1]*impulse;
    o1[2] += J1[2]*impulse;

    o2[0] -= J2[0]*impulse;
    o2[1] -= J2[1]*impulse;
    o2[2] -= J2[2]*impulse;
}

}

void Penalization::preventCollidingObstacles() const
{
  const auto & shapes = sim.obstacle_vector->getObstacleVector();
  const auto & infos  = sim.chiInfo();
  const size_t N = sim.obstacle_vector->nObstacles();
  sim.bCollisionID.clear();
  sim.bCollision = false;

  struct CollisionInfo // hitter and hittee, symmetry but we do things twice
  {
    Real iM = 0;
    Real iPosX = 0;
    Real iPosY = 0;
    Real iPosZ = 0;
    Real iMomX = 0;
    Real iMomY = 0;
    Real iMomZ = 0;
    Real ivecX = 0;
    Real ivecY = 0;
    Real ivecZ = 0;
    Real jM = 0;
    Real jPosX = 0;
    Real jPosY = 0;
    Real jPosZ = 0;
    Real jMomX = 0;
    Real jMomY = 0;
    Real jMomZ = 0;
    Real jvecX = 0;
    Real jvecY = 0;
    Real jvecZ = 0;
  };
  std::vector<CollisionInfo> collisions(N);

  std::vector <Real> n_vec(3*N,0.0);

  #pragma omp parallel for schedule(static)
  for (size_t i=0; i<N; ++i)
  for (size_t j=0; j<N; ++j)
  {
    if(i==j) continue;
    auto & coll = collisions[i];

    const auto& iBlocks = shapes[i]->obstacleBlocks;
    const Real iU0      = shapes[i]->transVel[0];
    const Real iU1      = shapes[i]->transVel[1];
    const Real iU2      = shapes[i]->transVel[2];
    const Real iomega0  = shapes[i]->angVel  [0];
    const Real iomega1  = shapes[i]->angVel  [1];
    const Real iomega2  = shapes[i]->angVel  [2];
    const Real iCx      = shapes[i]->centerOfMass[0];
    const Real iCy      = shapes[i]->centerOfMass[1];
    const Real iCz      = shapes[i]->centerOfMass[2];

    const auto& jBlocks = shapes[j]->obstacleBlocks;
    const Real jU0      = shapes[j]->transVel[0];
    const Real jU1      = shapes[j]->transVel[1];
    const Real jU2      = shapes[j]->transVel[2];
    const Real jomega0  = shapes[j]->angVel  [0];
    const Real jomega1  = shapes[j]->angVel  [1];
    const Real jomega2  = shapes[j]->angVel  [2];
    const Real jCx      = shapes[j]->centerOfMass[0];
    const Real jCy      = shapes[j]->centerOfMass[1];
    const Real jCz      = shapes[j]->centerOfMass[2];

    assert(iBlocks.size() == jBlocks.size());

    const size_t nBlocks = iBlocks.size();
    for (size_t k=0; k<nBlocks; ++k)
    {
      if ( iBlocks[k] == nullptr || jBlocks[k] == nullptr ) continue;

      const auto    & iSDF  = iBlocks[k]->sdfLab;
      const auto    & jSDF  = jBlocks[k]->sdfLab;
      const CHIMAT  & iChi  = iBlocks[k]->chi;
      const CHIMAT  & jChi  = jBlocks[k]->chi;
      const UDEFMAT & iUDEF = iBlocks[k]->udef;
      const UDEFMAT & jUDEF = jBlocks[k]->udef;

      for(int z=0; z<VectorBlock::sizeZ; ++z)
      for(int y=0; y<VectorBlock::sizeY; ++y)
      for(int x=0; x<VectorBlock::sizeX; ++x)
      {
        if(iChi[z][y][x] <= 0.0 || jChi[z][y][x] <= 0.0 ) continue;
        const auto pos = infos[k].pos<Real>(x, y, z);

        coll.iM    += iChi[z][y][x];
        coll.iPosX += iChi[z][y][x] * pos[0];
        coll.iPosY += iChi[z][y][x] * pos[1];
        coll.iPosZ += iChi[z][y][x] * pos[2];
        coll.iMomX += iChi[z][y][x] * (iU0 + iomega1* (pos[2] - iCz) - iomega2*(pos[1]-iCy) + iUDEF[z][y][x][0]);
        coll.iMomY += iChi[z][y][x] * (iU1 + iomega2* (pos[0] - iCx) - iomega0*(pos[2]-iCz) + iUDEF[z][y][x][1]);
        coll.iMomZ += iChi[z][y][x] * (iU2 + iomega0* (pos[1] - iCy) - iomega1*(pos[0]-iCx) + iUDEF[z][y][x][2]);
        coll.ivecX += iChi[z][y][x] * 0.5*(iSDF[z+1][y+1][x+2] - iSDF[z+1][y+1][x  ]);
        coll.ivecY += iChi[z][y][x] * 0.5*(iSDF[z+1][y+2][x+1] - iSDF[z+1][y  ][x+1]);
        coll.ivecZ += iChi[z][y][x] * 0.5*(iSDF[z+2][y+1][x+1] - iSDF[z  ][y+1][x+1]);

        coll.jM    += jChi[z][y][x];
        coll.jPosX += jChi[z][y][x] * pos[0];
        coll.jPosY += jChi[z][y][x] * pos[1];
        coll.jPosZ += jChi[z][y][x] * pos[2];
        coll.jMomX += jChi[z][y][x] * (jU0 + jomega1* (pos[2] - jCz) - jomega2*(pos[1]-jCy) + jUDEF[z][y][x][0]);
        coll.jMomY += jChi[z][y][x] * (jU1 + jomega2* (pos[0] - jCx) - jomega0*(pos[2]-jCz) + jUDEF[z][y][x][1]);
        coll.jMomZ += jChi[z][y][x] * (jU2 + jomega0* (pos[1] - jCy) - jomega1*(pos[0]-jCx) + jUDEF[z][y][x][2]);
        coll.jvecX += jChi[z][y][x] * 0.5*(jSDF[z+1][y+1][x+2] - jSDF[z+1][y+1][x  ]);
        coll.jvecY += jChi[z][y][x] * 0.5*(jSDF[z+1][y+2][x+1] - jSDF[z+1][y  ][x+1]);
        coll.jvecZ += jChi[z][y][x] * 0.5*(jSDF[z+2][y+1][x+1] - jSDF[z  ][y+1][x+1]);
      }
    }
  }

  std::vector<Real> buffer(20*N); //CollisionInfo holds 20 Reals
  for (size_t i = 0 ; i < N ; i++)
  {
    buffer[20*i     ] = collisions[i].iM   ;
    buffer[20*i + 1 ] = collisions[i].iPosX;
    buffer[20*i + 2 ] = collisions[i].iPosY;
    buffer[20*i + 3 ] = collisions[i].iPosZ;
    buffer[20*i + 4 ] = collisions[i].iMomX;
    buffer[20*i + 5 ] = collisions[i].iMomY;
    buffer[20*i + 6 ] = collisions[i].iMomZ;
    buffer[20*i + 7 ] = collisions[i].ivecX;
    buffer[20*i + 8 ] = collisions[i].ivecY;
    buffer[20*i + 9 ] = collisions[i].ivecZ;
    buffer[20*i + 10] = collisions[i].jM   ;
    buffer[20*i + 11] = collisions[i].jPosX;
    buffer[20*i + 12] = collisions[i].jPosY;
    buffer[20*i + 13] = collisions[i].jPosZ;
    buffer[20*i + 14] = collisions[i].jMomX;
    buffer[20*i + 15] = collisions[i].jMomY;
    buffer[20*i + 16] = collisions[i].jMomZ;
    buffer[20*i + 17] = collisions[i].jvecX;
    buffer[20*i + 18] = collisions[i].jvecY;
    buffer[20*i + 19] = collisions[i].jvecZ;
  }
  MPI_Allreduce(MPI_IN_PLACE, buffer.data(), buffer.size(), MPI_Real, MPI_SUM, sim.comm);
  for (size_t i = 0 ; i < N ; i++)
  {
    collisions[i].iM    = buffer[20*i     ];
    collisions[i].iPosX = buffer[20*i + 1 ];
    collisions[i].iPosY = buffer[20*i + 2 ];
    collisions[i].iPosZ = buffer[20*i + 3 ];
    collisions[i].iMomX = buffer[20*i + 4 ];
    collisions[i].iMomY = buffer[20*i + 5 ];
    collisions[i].iMomZ = buffer[20*i + 6 ];
    collisions[i].ivecX = buffer[20*i + 7 ];
    collisions[i].ivecY = buffer[20*i + 8 ];
    collisions[i].ivecZ = buffer[20*i + 9 ];
    collisions[i].jM    = buffer[20*i + 10];
    collisions[i].jPosX = buffer[20*i + 11];
    collisions[i].jPosY = buffer[20*i + 12];
    collisions[i].jPosZ = buffer[20*i + 13];
    collisions[i].jMomX = buffer[20*i + 14];
    collisions[i].jMomY = buffer[20*i + 15];
    collisions[i].jMomZ = buffer[20*i + 16];
    collisions[i].jvecX = buffer[20*i + 17];
    collisions[i].jvecY = buffer[20*i + 18];
    collisions[i].jvecZ = buffer[20*i + 19];
  }

  #pragma omp parallel for schedule(static)
  for (size_t i=0  ; i<N; ++i)
  for (size_t j=i+1; j<N; ++j)
  {
    const auto & coll       = collisions[i];
    const auto & coll_other = collisions[j];

    // less than one fluid element of overlap: wait to get closer. no hit
    if(coll.iM       < 8.0 || coll.jM       < 8.0) continue; //object i did not collide
    if(coll_other.iM < 8.0 || coll_other.jM < 8.0) continue; //object j did not collide

    if (std::fabs(coll.iPosX/coll.iM  - coll_other.iPosX/coll_other.iM ) > 0.2 ||
        std::fabs(coll.iPosY/coll.iM  - coll_other.iPosY/coll_other.iM ) > 0.2 ||
        std::fabs(coll.iPosZ/coll.iM  - coll_other.iPosZ/coll_other.iM ) > 0.2 ) //used 0.2 because fish lenght is 0.2 usually!
    {
        continue; // then both objects i and j collided, but not with each other!
    }

    // A collision happened!
    sim.bCollision = true;
    #pragma omp critical
	  {
	    sim.bCollisionID.push_back(i);
	    sim.bCollisionID.push_back(j);
	  }

    //1. Compute collision normal vector (NX,NY,NZ)
    const Real norm_i = std::sqrt(coll.ivecX*coll.ivecX + coll.ivecY*coll.ivecY + coll.ivecZ*coll.ivecZ);
    const Real norm_j = std::sqrt(coll.jvecX*coll.jvecX + coll.jvecY*coll.jvecY + coll.jvecZ*coll.jvecZ);
    const Real mX = coll.ivecX/norm_i - coll.jvecX/norm_j;
    const Real mY = coll.ivecY/norm_i - coll.jvecY/norm_j;
    const Real mZ = coll.ivecZ/norm_i - coll.jvecZ/norm_j;
    const Real inorm = 1.0/std::sqrt(mX*mX + mY*mY + mZ*mZ);
    const Real NX = mX * inorm;
    const Real NY = mY * inorm;
    const Real NZ = mZ * inorm;

    //If objects are already moving away from each other, don't do anything
    const Real hitVelX = coll.jMomX / coll.jM - coll.iMomX / coll.iM;
    const Real hitVelY = coll.jMomY / coll.jM - coll.iMomY / coll.iM;
    const Real hitVelZ = coll.jMomZ / coll.jM - coll.iMomZ / coll.iM;
    const Real projVel = hitVelX * NX + hitVelY * NY + hitVelZ * NZ;
    if(projVel<=0) continue; // vel goes away from collision: no need to bounce

    //2. Compute collision location
    const Real inv_iM = 1.0/coll.iM;
    const Real inv_jM = 1.0/coll.jM;
    const Real iPX = coll.iPosX * inv_iM; // object i collision location
    const Real iPY = coll.iPosY * inv_iM;
    const Real iPZ = coll.iPosZ * inv_iM;
    const Real jPX = coll.jPosX * inv_jM; // object j collision location
    const Real jPY = coll.jPosY * inv_jM;
    const Real jPZ = coll.jPosZ * inv_jM;
    const Real CX = 0.5*(iPX+jPX);
    const Real CY = 0.5*(iPY+jPY);
    const Real CZ = 0.5*(iPZ+jPZ);

    //3. Take care of the collision. Assume elastic collision (kinetic energy is conserved)
    const bool iforced = shapes[i]->bForcedInSimFrame[0] || shapes[i]->bForcedInSimFrame[1] || shapes[i]->bForcedInSimFrame[2];
    const bool jforced = shapes[j]->bForcedInSimFrame[0] || shapes[j]->bForcedInSimFrame[1] || shapes[j]->bForcedInSimFrame[2];
    const Real m1      = iforced ? 1e10*shapes[i]->mass : shapes[i]->mass;
    const Real m2      = jforced ? 1e10*shapes[j]->mass : shapes[j]->mass;
    //velocities of the two colliding points
    const Real vc1[3] = {coll.iMomX/coll.iM, coll.iMomY/coll.iM, coll.iMomZ/coll.iM};
    const Real vc2[3] = {coll.jMomX/coll.jM, coll.jMomY/coll.jM, coll.jMomZ/coll.jM};
    //save for printing
    const Real v1[3] = {shapes[i]->transVel[0],shapes[i]->transVel[1],shapes[i]->transVel[2]};
    const Real v2[3] = {shapes[j]->transVel[0],shapes[j]->transVel[1],shapes[j]->transVel[2]};
    const Real o1[3] = {shapes[i]->angVel  [0],shapes[i]->angVel  [1],shapes[i]->angVel  [2]};
    const Real o2[3] = {shapes[j]->angVel  [0],shapes[j]->angVel  [1],shapes[j]->angVel  [2]};

    ElasticCollision(m1,m2,&shapes[i]->J[0]           ,&shapes[j]->J[0]           ,
                           &shapes[i]->centerOfMass[0],&shapes[j]->centerOfMass[0],
                           NX,NY,NZ,CX,CY,CZ,vc1,vc2,
                           &shapes[i]->transVel[0],&shapes[j]->transVel[0],
                           &shapes[i]->angVel  [0],&shapes[j]->angVel  [0]);
    if ( sim.verbose )
    {
      #pragma omp critical
      {
        std::cout << "Collision between objects " << i << " and " << j << std::endl;
        std::cout << " iM   (0) = " << collisions[i].iM    << " jM   (1) = " << collisions[j].jM << std::endl;
        std::cout << " jM   (0) = " << collisions[i].jM    << " jM   (1) = " << collisions[j].iM << std::endl;
        std::cout << " Normal vector = (" << NX << "," << NY << "," << NZ << ")" << std::endl;
        std::cout << " Location      = (" << CX << "," << CY << "," << CZ << ")" << std::endl;
        std::cout << " Shape " << i << " before collision u    =(" <<  v1                 [0] << "," <<  v1                 [1] << "," <<  v1                 [2] << ")" << std::endl;
        std::cout << " Shape " << i << " after  collision u    =(" <<  shapes[i]->transVel[0] << "," <<  shapes[i]->transVel[1] << "," <<  shapes[i]->transVel[2] << ")" << std::endl;
        std::cout << " Shape " << j << " before collision u    =(" <<  v2                 [0] << "," <<  v2                 [1] << "," <<  v2                 [2] << ")" << std::endl;
        std::cout << " Shape " << j << " after  collision u    =(" <<  shapes[j]->transVel[0] << "," <<  shapes[j]->transVel[1] << "," <<  shapes[j]->transVel[2] << ")" << std::endl;
        std::cout << " Shape " << i << " before collision omega=(" <<  o1                 [0] << "," <<  o1                 [1] << "," <<  o1                 [2] << ")" << std::endl;
        std::cout << " Shape " << i << " after  collision omega=(" <<  shapes[i]->angVel  [0] << "," <<  shapes[i]->angVel  [1] << "," <<  shapes[i]->angVel  [2] << ")" << std::endl;
        std::cout << " Shape " << j << " before collision omega=(" <<  o2                 [0] << "," <<  o2                 [1] << "," <<  o2                 [2] << ")" << std::endl;
        std::cout << " Shape " << j << " after  collision omega=(" <<  shapes[j]->angVel  [0] << "," <<  shapes[j]->angVel  [1] << "," <<  shapes[j]->angVel  [2] << ")" << std::endl;
      }
    }
  }
}


Penalization::Penalization(SimulationData & s) : Operator(s) {}

void Penalization::operator()(const Real dt)
{
  if(sim.obstacle_vector->nObstacles() == 0) return;

  preventCollidingObstacles();

  std::vector<cubism::BlockInfo>& chiInfo = sim.chiInfo();
  std::vector<cubism::BlockInfo>& velInfo = sim.velInfo();
  #pragma omp parallel
  { // each thread needs to call its own non-const operator() function
    if(sim.bImplicitPenalization)
    {
      KernelPenalization<1> K(sim.dt, sim.lambda, sim.obstacle_vector);
      #pragma omp for schedule(dynamic, 1)
      for (size_t i = 0; i < chiInfo.size(); ++i) K(velInfo[i],chiInfo[i]);
    }
    else
    {
      KernelPenalization<0> K(sim.dt, sim.lambda, sim.obstacle_vector);
      #pragma omp for schedule(dynamic, 1)
      for (size_t i = 0; i < chiInfo.size(); ++i) K(velInfo[i],chiInfo[i]);
    }
  }

  kernelFinalizePenalizationForce(sim);
}

CubismUP_3D_NAMESPACE_END
