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

struct KernelPenalization
{
  const Real dt, invdt = 1.0/dt, lambda;
  const bool implicitPenalization;
  ObstacleVector * const obstacle_vector;

  KernelPenalization(const Real _dt, const Real _lambda, const bool _implicitPenalization, ObstacleVector* ov) :
    dt(_dt), lambda(_lambda), implicitPenalization(_implicitPenalization), obstacle_vector(ov) {}

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
      const Real X = implicitPenalization? (CHI[iz][iy][ix]>0.5?1.0:0.0) : CHI[iz][iy][ix];
      const Real penalFac = implicitPenalization? X*lambdaFac/(1+ X*lambdaFac*dt):X*lambdaFac;
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
    //Invert I
    const Real m00 = I[0];
    const Real m01 = I[3];
    const Real m02 = I[4];
    const Real m11 = I[1];
    const Real m12 = I[5];
    const Real m22 = I[2];
    Real a00 = m22*m11 - m12*m12;
    Real a01 = m02*m12 - m22*m01;
    Real a02 = m01*m12 - m02*m11;
    Real a11 = m22*m00 - m02*m02;
    Real a12 = m01*m02 - m00*m12;
    Real a22 = m00*m11 - m01*m01;
    const Real determinant =  1.0/((m00 * a00) + (m01 * a01) + (m02 * a02));
    a00 *= determinant;
    a01 *= determinant;
    a02 *= determinant;
    a11 *= determinant;
    a12 *= determinant;
    a22 *= determinant;
    const Real aux_0 = ( Rc[1] - R[1] )*N[2] - ( Rc[2] - R[2] )*N[1];
    const Real aux_1 = ( Rc[2] - R[2] )*N[0] - ( Rc[0] - R[0] )*N[2];
    const Real aux_2 = ( Rc[0] - R[0] )*N[1] - ( Rc[1] - R[1] )*N[0];
    J[0] = a00*aux_0 + a01*aux_1 + a02*aux_2;
    J[1] = a01*aux_0 + a11*aux_1 + a12*aux_2;
    J[2] = a02*aux_0 + a12*aux_1 + a22*aux_2;
}

void ElasticCollision(const Real  m1,const Real  m2,const Real *I1,const Real *I2,
                      const Real *v1,const Real *v2,const Real *o1,const Real *o2,
                      const Real *C1,const Real *C2,const Real  NX,const Real  NY,const Real  NZ,
                      const Real  CX,const Real  CY,const Real  CZ,const Real *vc1,const Real *vc2,
                      Real *hv1,Real *hv2,Real *ho1,Real *ho2)
{
    const Real e = 1.0; // coefficient of restitution
    const Real N[3] ={NX,NY,NZ};
    const Real C[3] ={CX,CY,CZ};
    const Real k1[3] = { N[0]/m1, N[1]/m1, N[2]/m1};
    const Real k2[3] = {-N[0]/m2,-N[1]/m2,-N[2]/m2};
    Real J1[3];
    Real J2[3]; 
    ComputeJ(C,C1,N,I1,J1);
    ComputeJ(C,C2,N,I2,J2);
    const Real nom = (e+1)*((vc1[0]-vc2[0])*N[0] + (vc1[1]-vc2[1])*N[1] + (vc1[2]-vc2[2])*N[2]);
    const Real denom = -(1.0/m1+1.0/m2) + 
               -( ( J1[1]*(C[2]-C1[2]) - J1[2]*(C[1]-C1[1]) ) *N[0]+
                  ( J1[2]*(C[0]-C1[0]) - J1[0]*(C[2]-C1[2]) ) *N[1]+
                  ( J1[0]*(C[1]-C1[1]) - J1[1]*(C[0]-C1[0]) ) *N[2])
               -( ( J2[1]*(C[2]-C2[2]) - J2[2]*(C[1]-C2[1]) ) *N[0]+
                  ( J2[2]*(C[0]-C2[0]) - J2[0]*(C[2]-C2[2]) ) *N[1]+
                  ( J2[0]*(C[1]-C2[1]) - J2[1]*(C[0]-C2[0]) ) *N[2]);
    const Real impulse = nom/(denom+1e-21);
    hv1[0] = v1[0] + k1[0]*impulse;
    hv1[1] = v1[1] + k1[1]*impulse;
    hv1[2] = v1[2] + k1[2]*impulse;
    hv2[0] = v2[0] + k2[0]*impulse;
    hv2[1] = v2[1] + k2[1]*impulse;
    hv2[2] = v2[2] + k2[2]*impulse;
    ho1[0] = o1[0] + J1[0]*impulse;
    ho1[1] = o1[1] + J1[1]*impulse;
    ho1[2] = o1[2] + J1[2]*impulse;
    ho2[0] = o2[0] - J2[0]*impulse;
    ho2[1] = o2[1] - J2[1]*impulse;
    ho2[2] = o2[2] - J2[2]*impulse;
}

}

void Penalization::preventCollidingObstacles() const
{
    const auto & shapes = sim.obstacle_vector->getObstacleVector();
    const auto & infos  = sim.chiInfo();
    const size_t N = sim.obstacle_vector->nObstacles();
    sim.bCollisionID.clear();

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

    #pragma omp parallel for schedule(static)
    for (size_t i=0; i<N; ++i)
    {
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

        for (size_t j=0; j<N; ++j)
        {
            if(i==j) continue;
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

            Real imagmax = 0.0;
            Real jmagmax = 0.0;

            assert(iBlocks.size() == jBlocks.size());
            for (size_t k=0; k<iBlocks.size(); ++k)
            {
                if ( iBlocks[k] == nullptr || jBlocks[k] == nullptr ) continue;

                const auto & iSDF  = iBlocks[k]->sdfLab;
                const auto & jSDF  = jBlocks[k]->sdfLab;
                const auto & iChi  = iBlocks[k]->chi;
                const auto & jChi  = jBlocks[k]->chi;
                const auto & iUDEF = iBlocks[k]->udef;
                const auto & jUDEF = jBlocks[k]->udef;

                for(int z=0; z<VectorBlock::sizeZ; ++z)
                for(int y=0; y<VectorBlock::sizeY; ++y)
                for(int x=0; x<VectorBlock::sizeX; ++x)
                {
                    if(iChi[z][y][x] <= 0.0 || jChi[z][y][x] <= 0.0) continue;

                    const auto p = infos[k].pos<Real>(x,y,z);
                    const Real iMomX = iU0 + iomega1*(p[2]-iCz) - iomega2*(p[1]-iCy) + iUDEF[z][y][x][0];
                    const Real iMomY = iU1 + iomega2*(p[0]-iCx) - iomega0*(p[2]-iCz) + iUDEF[z][y][x][1];
                    const Real iMomZ = iU2 + iomega0*(p[1]-iCy) - iomega1*(p[0]-iCx) + iUDEF[z][y][x][2];
                    const Real jMomX = jU0 + jomega1*(p[2]-jCz) - jomega2*(p[1]-jCy) + jUDEF[z][y][x][0];
                    const Real jMomY = jU1 + jomega2*(p[0]-jCx) - jomega0*(p[2]-jCz) + jUDEF[z][y][x][1];
                    const Real jMomZ = jU2 + jomega0*(p[1]-jCy) - jomega1*(p[0]-jCx) + jUDEF[z][y][x][2];

                    const Real imag = iMomX*iMomX+iMomY*iMomY+iMomZ*iMomZ;
                    const Real jmag = jMomX*jMomX+jMomY*jMomY+jMomZ*jMomZ;

                    const Real ivecX = iSDF[z+1][y+1][x+2] - iSDF[z+1][y+1][x  ];
                    const Real ivecY = iSDF[z+1][y+2][x+1] - iSDF[z+1][y  ][x+1];
                    const Real ivecZ = iSDF[z+2][y+1][x+1] - iSDF[z  ][y+1][x+1];
                    const Real jvecX = jSDF[z+1][y+1][x+2] - jSDF[z+1][y+1][x  ];
                    const Real jvecY = jSDF[z+1][y+2][x+1] - jSDF[z+1][y  ][x+1];
                    const Real jvecZ = jSDF[z+2][y+1][x+1] - jSDF[z  ][y+1][x+1];
                    const Real normi = 1.0/(sqrt(ivecX*ivecX+ivecY*ivecY+ivecZ*ivecZ)+1e-21);
                    const Real normj = 1.0/(sqrt(jvecX*jvecX+jvecY*jvecY+jvecZ*jvecZ)+1e-21);

                    coll.iM    += 1;
                    coll.iPosX += p[0];
                    coll.iPosY += p[1];
                    coll.iPosZ += p[2];
                    coll.ivecX += ivecX*normi;
                    coll.ivecY += ivecY*normi;
                    coll.ivecZ += ivecZ*normi;
                    if (imag > imagmax)
                    {
                        imagmax = imag;
                        coll.iMomX = iMomX;
                        coll.iMomY = iMomY;
                        coll.iMomZ = iMomZ;
                    }

                    coll.jM    += 1;
                    coll.jPosX += p[0];
                    coll.jPosY += p[1];
                    coll.jPosZ += p[2];
                    coll.jvecX += jvecX*normj;
                    coll.jvecY += jvecY*normj;
                    coll.jvecZ += jvecZ*normj;
                    if (jmag > jmagmax)
                    {
                        jmagmax = jmag;
                        coll.jMomX = jMomX;
                        coll.jMomY = jMomY;
                        coll.jMomZ = jMomZ;
                    }
                }
            }
        }
    }

    std::vector<Real> buffer(20*N); //CollisionInfo holds 20 Reals
    std::vector<Real> buffermax(2*N);
    for (size_t i = 0 ; i < N ; i++)
    {
        const auto & coll = collisions[i];
        buffermax[2*i  ] = coll.iMomX*coll.iMomX+ coll.iMomY*coll.iMomY+coll.iMomZ*coll.iMomZ;
        buffermax[2*i+1] = coll.jMomX*coll.jMomX+ coll.jMomY*coll.jMomY+coll.jMomZ*coll.jMomZ;
    }
    MPI_Allreduce(MPI_IN_PLACE, buffermax.data(), buffermax.size(), MPI_Real, MPI_MAX, sim.comm);

    for (size_t i = 0 ; i < N ; i++)
    {
        const auto & coll = collisions[i];
        const Real maxi = coll.iMomX*coll.iMomX+ coll.iMomY*coll.iMomY+coll.iMomZ*coll.iMomZ;
        const Real maxj = coll.jMomX*coll.jMomX+ coll.jMomY*coll.jMomY+coll.jMomZ*coll.jMomZ;
        const bool iok = std::fabs(maxi - buffermax[2*i  ]) < 1e-10;
        const bool jok = std::fabs(maxj - buffermax[2*i+1]) < 1e-10;
        buffer[20*i     ] = coll.iM   ;
        buffer[20*i + 1 ] = coll.iPosX;
        buffer[20*i + 2 ] = coll.iPosY;
        buffer[20*i + 3 ] = coll.iPosZ;
        buffer[20*i + 4 ] = iok ? coll.iMomX : 0;
        buffer[20*i + 5 ] = iok ? coll.iMomY : 0;
        buffer[20*i + 6 ] = iok ? coll.iMomZ : 0;
        buffer[20*i + 7 ] = coll.ivecX;
        buffer[20*i + 8 ] = coll.ivecY;
        buffer[20*i + 9 ] = coll.ivecZ;
        buffer[20*i + 10] = coll.jM   ;
        buffer[20*i + 11] = coll.jPosX;
        buffer[20*i + 12] = coll.jPosY;
        buffer[20*i + 13] = coll.jPosZ;
        buffer[20*i + 14] = jok ? coll.jMomX :0;
        buffer[20*i + 15] = jok ? coll.jMomY :0;
        buffer[20*i + 16] = jok ? coll.jMomZ :0;
        buffer[20*i + 17] = coll.jvecX;
        buffer[20*i + 18] = coll.jvecY;
        buffer[20*i + 19] = coll.jvecZ;
    }
    MPI_Allreduce(MPI_IN_PLACE, buffer.data(), buffer.size(), MPI_Real, MPI_SUM, sim.comm);
    for (size_t i = 0 ; i < N ; i++)
    {
        auto & coll = collisions[i];
        coll.iM    = buffer[20*i     ];
        coll.iPosX = buffer[20*i + 1 ];
        coll.iPosY = buffer[20*i + 2 ];
        coll.iPosZ = buffer[20*i + 3 ];
        coll.iMomX = buffer[20*i + 4 ];
        coll.iMomY = buffer[20*i + 5 ];
        coll.iMomZ = buffer[20*i + 6 ];
        coll.ivecX = buffer[20*i + 7 ];
        coll.ivecY = buffer[20*i + 8 ];
        coll.ivecZ = buffer[20*i + 9 ];
        coll.jM    = buffer[20*i + 10];
        coll.jPosX = buffer[20*i + 11];
        coll.jPosY = buffer[20*i + 12];
        coll.jPosZ = buffer[20*i + 13];
        coll.jMomX = buffer[20*i + 14];
        coll.jMomY = buffer[20*i + 15];
        coll.jMomZ = buffer[20*i + 16];
        coll.jvecX = buffer[20*i + 17];
        coll.jvecY = buffer[20*i + 18];
        coll.jvecZ = buffer[20*i + 19];
    }

    #pragma omp parallel for schedule(static)
    for (size_t i=0  ; i<N; ++i)
    for (size_t j=i+1; j<N; ++j)
    {
        const Real m1 = shapes[i]->mass;
        const Real m2 = shapes[j]->mass;
        const Real v1[3]={shapes[i]->transVel[0],shapes[i]->transVel[1],shapes[i]->transVel[2]};
        const Real o1[3]={shapes[i]->  angVel[0],shapes[i]->  angVel[1],shapes[i]->  angVel[2]};
        const Real v2[3]={shapes[j]->transVel[0],shapes[j]->transVel[1],shapes[j]->transVel[2]};
        const Real o2[3]={shapes[j]->  angVel[0],shapes[j]->  angVel[1],shapes[j]->  angVel[2]};
        const Real I1[6]={shapes[i]->J[0],shapes[i]->J[1],shapes[i]->J[2],shapes[i]->J[3],shapes[i]->J[4],shapes[i]->J[5]};
        const Real I2[6]={shapes[j]->J[0],shapes[j]->J[1],shapes[j]->J[2],shapes[j]->J[3],shapes[j]->J[4],shapes[j]->J[5]};
        const Real C1[3]={shapes[i]->centerOfMass[0],shapes[i]->centerOfMass[1],shapes[i]->centerOfMass[2]};
        const Real C2[3]={shapes[j]->centerOfMass[0],shapes[j]->centerOfMass[1],shapes[j]->centerOfMass[2]};

        auto & coll       = collisions[i];
        auto & coll_other = collisions[j];

        // less than 'tolerance' fluid element(s) of overlap: wait to get closer. no hit
        const Real tolerance = 0.001;
        if(coll.iM       < tolerance || coll.jM       < tolerance) continue;
        if(coll_other.iM < tolerance || coll_other.jM < tolerance) continue;

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
        const Real projVel = (coll.jMomX  - coll.iMomX) * NX + (coll.jMomY  - coll.iMomY) * NY + (coll.jMomZ  - coll.iMomZ) * NZ;
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
        const Real vc1[3] = {coll.iMomX, coll.iMomY, coll.iMomZ};
        const Real vc2[3] = {coll.jMomX, coll.jMomY, coll.jMomZ};
        Real ho1[3];
        Real ho2[3];
        Real hv1[3];
        Real hv2[3];
        const bool iforced = shapes[i]->bForcedInSimFrame[0] || shapes[i]->bForcedInSimFrame[1] || shapes[i]->bForcedInSimFrame[2];
        const bool jforced = shapes[j]->bForcedInSimFrame[0] || shapes[j]->bForcedInSimFrame[1] || shapes[j]->bForcedInSimFrame[2];
        const Real m1_i = iforced ? 1e10*m1 : m1;
        const Real m2_j = jforced ? 1e10*m2 : m2;
        ElasticCollision(m1_i,m2_j,I1,I2,v1,v2,o1,o2,C1,C2,NX,NY,NZ,CX,CY,CZ,vc1,vc2,hv1,hv2,ho1,ho2);

        shapes[i]->transVel[0] = hv1[0];
        shapes[i]->transVel[1] = hv1[1];
        shapes[i]->transVel[2] = hv1[2];
        shapes[j]->transVel[0] = hv2[0];
        shapes[j]->transVel[1] = hv2[1];
        shapes[j]->transVel[2] = hv2[2];
        shapes[i]->angVel[0] = ho1[0];
        shapes[i]->angVel[1] = ho1[1];
        shapes[i]->angVel[2] = ho1[2];
        shapes[j]->angVel[0] = ho2[0];
        shapes[j]->angVel[1] = ho2[1];
        shapes[j]->angVel[2] = ho2[2];

	    shapes[i]-> u_collision = hv1[0];
        shapes[i]-> v_collision = hv1[1];
        shapes[i]-> w_collision = hv1[2];
        shapes[i]->ox_collision = ho1[0];
        shapes[i]->oy_collision = ho1[1];
        shapes[i]->oz_collision = ho1[2];
        shapes[j]-> u_collision = hv2[0];
        shapes[j]-> v_collision = hv2[1];
        shapes[j]-> w_collision = hv2[2];
        shapes[j]->ox_collision = ho2[0];
        shapes[j]->oy_collision = ho2[1];
        shapes[j]->oz_collision = ho2[2];
        shapes[i]->collision_counter = 0.01*sim.dt;
        shapes[j]->collision_counter = 0.01*sim.dt;

        if ( sim.verbose )
        {
            #pragma omp critical
            {
                std::cout << "Collision between objects " << i << " and " << j << std::endl;
                std::cout << " iM   (0) = " << collisions[i].iM    << " jM   (1) = " << collisions[j].jM << std::endl;
                std::cout << " jM   (0) = " << collisions[i].jM    << " jM   (1) = " << collisions[j].iM << std::endl;
                std::cout << " Normal vector = (" << NX << "," << NY << "," << NZ << ")" << std::endl;
                std::cout << " Location      = (" << CX << "," << CY << "," << CZ << ")" << std::endl;
                std::cout << " Shape " << i << " before collision u    =(" <<  v1[0] << "," <<  v1[1] << "," <<  v1[2] << ")" << std::endl;
                std::cout << " Shape " << i << " after  collision u    =(" << hv1[0] << "," << hv1[1] << "," << hv1[2] << ")" << std::endl;
                std::cout << " Shape " << j << " before collision u    =(" <<  v2[0] << "," <<  v2[1] << "," <<  v2[2] << ")" << std::endl;
                std::cout << " Shape " << j << " after  collision u    =(" << hv2[0] << "," << hv2[1] << "," << hv2[2] << ")" << std::endl;
                std::cout << " Shape " << i << " before collision omega=(" <<  o1[0] << "," <<  o1[1] << "," <<  o1[2] << ")" << std::endl;
                std::cout << " Shape " << i << " after  collision omega=(" << ho1[0] << "," << ho1[1] << "," << ho1[2] << ")" << std::endl;
                std::cout << " Shape " << j << " before collision omega=(" <<  o2[0] << "," <<  o2[1] << "," <<  o2[2] << ")" << std::endl;
                std::cout << " Shape " << j << " after  collision omega=(" << ho2[0] << "," << ho2[1] << "," << ho2[2] << ")" << std::endl;
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
  {
    // each thread needs to call its own non-const operator() function
    KernelPenalization K(dt, sim.lambda, sim.bImplicitPenalization, sim.obstacle_vector);
    #pragma omp for schedule(dynamic, 1)
    for (size_t i = 0; i < chiInfo.size(); ++i) K(velInfo[i],chiInfo[i]);
  }

  kernelFinalizePenalizationForce(sim);
}

CubismUP_3D_NAMESPACE_END
