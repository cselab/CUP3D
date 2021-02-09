//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#include "Penalization.h"
#include "../obstacles/ObstacleVector.h"

CubismUP_3D_NAMESPACE_BEGIN
using namespace cubism;

namespace {

using CHIMAT = Real[CUP_BLOCK_SIZE][CUP_BLOCK_SIZE][CUP_BLOCK_SIZE];
using UDEFMAT = Real[CUP_BLOCK_SIZE][CUP_BLOCK_SIZE][CUP_BLOCK_SIZE][3];

template<bool implicitPenalization>
struct KernelPenalization : public ObstacleVisitor
{
  const Real dt, invdt = 1.0/dt, lambda;
  ObstacleVector * const obstacle_vector;
  const cubism::BlockInfo * info_ptr = nullptr;

  KernelPenalization(double _dt, double _lambda, ObstacleVector* ov) :
    dt(_dt), lambda(_lambda), obstacle_vector(ov) {}

  void operator()(const cubism::BlockInfo& info)
  {
    // first store the lab and info, then do visitor
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

    const CHIMAT & __restrict__ CHI = o->chi;
    const UDEFMAT & __restrict__ UDEF = o->udef;
    FluidBlock& b = *(FluidBlock*)info.ptrBlock;
    const std::array<double,3> CM = obstacle->getCenterOfMass();
    const std::array<double,3> vel = obstacle->getTranslationVelocity();
    const std::array<double,3> omega = obstacle->getAngularVelocity();
    const Real dv = std::pow(info.h_gridpoint, 3);

    // Obstacle-specific lambda, useful for gradually adding an obstacle to the flow.
    const double rampUp = obstacle->lambda_factor;
    // lambda = 1/dt hardcoded for expl time int, other options are wrong.
    const double lambdaFac = rampUp * (implicitPenalization? lambda : invdt);

    double &FX = o->FX, &FY = o->FY, &FZ = o->FZ;
    double &TX = o->TX, &TY = o->TY, &TZ = o->TZ;
    FX = 0; FY = 0; FZ = 0; TX = 0; TY = 0; TZ = 0;

    for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
    for(int iy=0; iy<FluidBlock::sizeY; ++iy)
    for(int ix=0; ix<FluidBlock::sizeX; ++ix)
    {
      // What if multiple obstacles share a block? Do not write udef onto
      // grid if CHI stored on the grid is greater than obst's CHI.
      if(b(ix,iy,iz).chi > CHI[iz][iy][ix]) continue;
      if(CHI[iz][iy][ix] <= 0) continue; // no need to do anything
      double p[3]; info.pos(p, ix, iy, iz);
      p[0] -= CM[0]; p[1] -= CM[1]; p[2] -= CM[2];

      const double U_TOT[3] = {
          vel[0] + omega[1]*p[2] - omega[2]*p[1] + UDEF[iz][iy][ix][0],
          vel[1] + omega[2]*p[0] - omega[0]*p[2] + UDEF[iz][iy][ix][1],
          vel[2] + omega[0]*p[1] - omega[1]*p[0] + UDEF[iz][iy][ix][2]
      };
      //chi is either 0 or 1 (not approximated by a mollified Heaviside here
      //as this is a pointwise discrete scheme
      const Real X = CHI[iz][iy][ix];
      const Real penalFac = implicitPenalization? X*lambdaFac/(1+X*lambdaFac*dt):X*lambdaFac;

      const Real FPX = penalFac * (U_TOT[0] - b(ix,iy,iz).u);
      const Real FPY = penalFac * (U_TOT[1] - b(ix,iy,iz).v);
      const Real FPZ = penalFac * (U_TOT[2] - b(ix,iy,iz).w);
      // What if two obstacles overlap? Let's plus equal. We will need a
      // repulsion term of the velocity at some point in the code.
      b(ix,iy,iz).u = b(ix,iy,iz).u + dt * FPX;
      b(ix,iy,iz).v = b(ix,iy,iz).v + dt * FPY;
      b(ix,iy,iz).w = b(ix,iy,iz).w + dt * FPZ;

      FX += dv * FPX; FY += dv * FPY; FZ += dv * FPZ;
      TX += dv * ( p[1] * FPZ - p[2] * FPY );
      TY += dv * ( p[2] * FPX - p[0] * FPZ );
      TZ += dv * ( p[0] * FPY - p[1] * FPX );
    }
  }
};

struct KernelFinalizePenalizationForce : public ObstacleVisitor
{
  FluidGridMPI * const grid;

  KernelFinalizePenalizationForce(FluidGridMPI*g) : grid(g) { }

  void visit(Obstacle* const obst)
  {
    static constexpr int nQoI = 6;
    double M[nQoI] = { 0 };
    const auto& oBlock = obst->getObstacleBlocks();
    #pragma omp parallel for schedule(static) reduction(+ : M[:nQoI])
    for (size_t i=0; i<oBlock.size(); ++i) {
      if(oBlock[i] == nullptr) continue;
      M[0] += oBlock[i]->FX; M[1] += oBlock[i]->FY; M[2] += oBlock[i]->FZ;
      M[3] += oBlock[i]->TX; M[4] += oBlock[i]->TY; M[5] += oBlock[i]->TZ;
    }
    const auto comm = grid->getCartComm();
    MPI_Allreduce(MPI_IN_PLACE, M, nQoI, MPI_DOUBLE, MPI_SUM, comm);
    obst->force[0]  = M[0]; obst->force[1]  = M[1]; obst->force[2]  = M[2];
    obst->torque[0] = M[3]; obst->torque[1] = M[4]; obst->torque[2] = M[5];
  }
};

}

void Penalization::preventCollidingObstacles() const
{
    using CHI_MAT = Real[FluidBlock::sizeZ][FluidBlock::sizeY][FluidBlock::sizeX];
    using UDEFMAT = Real[FluidBlock::sizeZ][FluidBlock::sizeY][FluidBlock::sizeX][3];
    const Real EPS = 1e-50;

    const auto & shapes = sim.obstacle_vector->getObstacleVector();
    const auto & infos  = sim.grid->getBlocksInfo();
    const size_t N = sim.obstacle_vector->nObstacles();

    //if (shapes[0]->bForcedInSimFrame[0]) shapes[0]->bForcedInSimFrame[0] = false;
    //if (shapes[0]->bForcedInSimFrame[1]) shapes[0]->bForcedInSimFrame[1] = false;
    //if (shapes[0]->bForcedInSimFrame[2]) shapes[0]->bForcedInSimFrame[2] = false;
    //if (shapes[0]->bFixFrameOfRef   [0]) shapes[0]->bFixFrameOfRef   [0] = false;
    //if (shapes[0]->bFixFrameOfRef   [1]) shapes[0]->bFixFrameOfRef   [1] = false;
    //if (shapes[0]->bFixFrameOfRef   [2]) shapes[0]->bFixFrameOfRef   [2] = false;

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

    std::vector <double> n_vec(3*N,0.0);

    #pragma omp parallel for schedule(static)
    for (size_t i=0; i<N; ++i)
    for (size_t j=0; j<N; ++j)
    {
        if(i==j) continue;
        auto & coll = collisions[i];

        const auto& iBlocks = shapes[i]->obstacleBlocks;
        const Real iU0      = shapes[i]->transVel[0];
        const Real iU1      = shapes[i]->transVel[1];
        const Real iU2      = shapes[i]->transVel[1];
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

            const CHI_MAT & iSDF  = iBlocks[k]->sdf;
            const CHI_MAT & jSDF  = jBlocks[k]->sdf;

            const CHI_MAT & iChi  = iBlocks[k]->chi;
            const CHI_MAT & jChi  = jBlocks[k]->chi;

            const UDEFMAT & iUDEF = iBlocks[k]->udef;
            const UDEFMAT & jUDEF = jBlocks[k]->udef;

            for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
            for(int iy=0; iy<FluidBlock::sizeY; ++iy)
            for(int ix=0; ix<FluidBlock::sizeX; ++ix)
            {

                if(iChi[iz][iy][ix] <= 0.0 || jChi[iz][iy][ix] <= 0.0 ) continue;

                const auto pos = infos[k].pos<Real>(ix, iy, iz);

                const Real iUr0 = iomega1* (pos[2] - iCz) - iomega2*(pos[1]-iCy);
                const Real iUr1 = iomega2* (pos[0] - iCx) - iomega0*(pos[2]-iCz);
                const Real iUr2 = iomega0* (pos[1] - iCy) - iomega1*(pos[0]-iCx);
                coll.iM    += iChi[iz][iy][ix];
                coll.iPosX += iChi[iz][iy][ix] * pos[0];
                coll.iPosY += iChi[iz][iy][ix] * pos[1];
                coll.iPosZ += iChi[iz][iy][ix] * pos[2];
                coll.iMomX += iChi[iz][iy][ix] * (iU0 + iUr0 + iUDEF[iz][iy][ix][0]);
                coll.iMomY += iChi[iz][iy][ix] * (iU1 + iUr1 + iUDEF[iz][iy][ix][1]);
                coll.iMomZ += iChi[iz][iy][ix] * (iU2 + iUr2 + iUDEF[iz][iy][ix][2]);

                const Real jUr0 = jomega1* (pos[2] - jCz) - iomega2*(pos[1]-jCy);
                const Real jUr1 = jomega2* (pos[0] - jCx) - iomega0*(pos[2]-jCz);
                const Real jUr2 = jomega0* (pos[1] - jCy) - iomega1*(pos[0]-jCx);
                coll.jM    += jChi[iz][iy][ix];
                coll.jPosX += jChi[iz][iy][ix] * pos[0];
                coll.jPosY += jChi[iz][iy][ix] * pos[1];
                coll.jPosZ += jChi[iz][iy][ix] * pos[2];
                coll.jMomX += jChi[iz][iy][ix] * (jU0 + jUr0 + jUDEF[iz][iy][ix][0]);
                coll.jMomY += jChi[iz][iy][ix] * (jU1 + jUr1 + jUDEF[iz][iy][ix][1]);
                coll.jMomZ += jChi[iz][iy][ix] * (jU2 + jUr2 + jUDEF[iz][iy][ix][1]);

                coll.ivecX += //iChi[iz][iy][ix] *
                 ( (ix == 0) ? (iSDF[iz][iy][ix+1] - iSDF[iz][iy][ix]) : ( (ix == FluidBlock::sizeX-1) ? (iSDF[iz][iy][ix] - iSDF[iz][iy][ix-1]) : 0.5*(iSDF[iz][iy][ix+1] - iSDF[iz][iy][ix-1]) ) );
                coll.ivecY += //iChi[iz][iy][ix] *
                 ( (iy == 0) ? (iSDF[iz][iy+1][ix] - iSDF[iz][iy][ix]) : ( (iy == FluidBlock::sizeY-1) ? (iSDF[iz][iy][ix] - iSDF[iz][iy-1][ix]) : 0.5*(iSDF[iz][iy+1][ix] - iSDF[iz][iy-1][ix]) ) );
                coll.ivecZ += //iChi[iz][iy][ix] *
                 ( (iz == 0) ? (iSDF[iz+1][iy][ix] - iSDF[iz][iy][ix]) : ( (iz == FluidBlock::sizeZ-1) ? (iSDF[iz][iy][ix] - iSDF[iz-1][iy][ix]) : 0.5*(iSDF[iz+1][iy][ix] - iSDF[iz-1][iy][ix]) ) );

                coll.jvecX += //jChi[iz][iy][ix] *
                 ( (ix == 0) ? (jSDF[iz][iy][ix+1] - jSDF[iz][iy][ix]) : ( (ix == FluidBlock::sizeX-1) ? (jSDF[iz][iy][ix] - jSDF[iz][iy][ix-1]) : 0.5*(jSDF[iz][iy][ix+1] - jSDF[iz][iy][ix-1]) ) );
                coll.jvecY += //jChi[iz][iy][ix] *
                 ( (iy == 0) ? (jSDF[iz][iy+1][ix] - jSDF[iz][iy][ix]) : ( (iy == FluidBlock::sizeY-1) ? (jSDF[iz][iy][ix] - jSDF[iz][iy-1][ix]) : 0.5*(jSDF[iz][iy+1][ix] - jSDF[iz][iy-1][ix]) ) );
                coll.jvecZ += //jChi[iz][iy][ix] *
                 ( (iz == 0) ? (jSDF[iz+1][iy][ix] - jSDF[iz][iy][ix]) : ( (iz == FluidBlock::sizeZ-1) ? (jSDF[iz][iy][ix] - jSDF[iz-1][iy][ix]) : 0.5*(jSDF[iz+1][iy][ix] - jSDF[iz-1][iy][ix]) ) );
            }
        }
    }

    std::vector<double> buffer(20*N); //CollisionInfo holds 20 doubles
    for (size_t i = 0 ; i < N ; i++)
    {
        auto & coll = collisions[i];
        buffer[20*i     ] = coll.iM   ;
        buffer[20*i + 1 ] = coll.iPosX;
        buffer[20*i + 2 ] = coll.iPosY;
        buffer[20*i + 3 ] = coll.iPosZ;
        buffer[20*i + 4 ] = coll.iMomX;
        buffer[20*i + 5 ] = coll.iMomY;
        buffer[20*i + 6 ] = coll.iMomZ;
        buffer[20*i + 7 ] = coll.ivecX;
        buffer[20*i + 8 ] = coll.ivecY;
        buffer[20*i + 9 ] = coll.ivecZ;
        buffer[20*i + 10] = coll.jM   ;
        buffer[20*i + 11] = coll.jPosX;
        buffer[20*i + 12] = coll.jPosY;
        buffer[20*i + 13] = coll.jPosZ;
        buffer[20*i + 14] = coll.jMomX;
        buffer[20*i + 15] = coll.jMomY;
        buffer[20*i + 16] = coll.jMomZ;
        buffer[20*i + 17] = coll.jvecX;
        buffer[20*i + 18] = coll.jvecY;
        buffer[20*i + 19] = coll.jvecZ;

    }
    MPI_Allreduce(MPI_IN_PLACE, buffer.data(), buffer.size(), MPI_DOUBLE, MPI_SUM, grid->getCartComm());
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
    for (size_t i=0; i<N; ++i)
    for (size_t j=0; j<N; ++j)
    {
        if (i==j) continue;

        auto & coll = collisions[i];
        // less than one fluid element of overlap: wait to get closer. no hit
        if(coll.iM < 20 || coll.jM < 20) continue;

        //1. Compute collision normal vector (NX,NY,NZ)
        const Real inv_iM = 1.0/coll.iM;
        const Real inv_jM = 1.0/coll.jM;
        const Real  inorm = 1.0/std::sqrt(coll.ivecX*coll.ivecX+coll.ivecY*coll.ivecY+coll.ivecZ*coll.ivecZ);
        const Real NX = coll.ivecX * inorm;
        const Real NY = coll.ivecY * inorm;
        const Real NZ = coll.ivecZ * inorm;

        //2. Compute collision location
        const Real iPX = coll.iPosX * inv_iM; // object i collision location
        const Real iPY = coll.iPosY * inv_iM;
        const Real iPZ = coll.iPosZ * inv_iM;
        const Real jPX = coll.jPosX * inv_jM; // object j collision location
        const Real jPY = coll.jPosY * inv_jM;
        const Real jPZ = coll.jPosZ * inv_jM;
        const Real CX = 0.5*(iPX+jPX);
        const Real CY = 0.5*(iPY+jPY);
        const Real CZ = 0.5*(iPZ+jPZ);

        //3. Compute collision velocity (hitVelX,hitVelY,hitVelZ)
        const Real hitVelX = coll.jMomX / coll.jM - coll.iMomX / coll.iM;
        const Real hitVelY = coll.jMomY / coll.jM - coll.iMomY / coll.iM;
        const Real hitVelZ = coll.jMomZ / coll.jM - coll.iMomZ / coll.iM;

        //4. Project velocity to collision normal direction
        const Real projVel = hitVelX * NX + hitVelY * NY + hitVelZ * NZ;

        if(projVel<=0) continue; // vel goes away from collision: no need to bounce

        std::cout << "Collision between objects " << i << " and " << j << std::endl;
        std::cout << " collision velocity = " << projVel << std::endl;

        //if (shapes[i]->bForcedInSimFrame[0]) shapes[i]->bForcedInSimFrame[0] = false;
        //if (shapes[i]->bForcedInSimFrame[1]) shapes[i]->bForcedInSimFrame[1] = false;
        //if (shapes[i]->bForcedInSimFrame[2]) shapes[i]->bForcedInSimFrame[2] = false;
        //if (shapes[j]->bForcedInSimFrame[0]) shapes[j]->bForcedInSimFrame[0] = false;
        //if (shapes[j]->bForcedInSimFrame[1]) shapes[j]->bForcedInSimFrame[1] = false;
        //if (shapes[j]->bForcedInSimFrame[2]) shapes[j]->bForcedInSimFrame[2] = false;

        //5. Take care of the collision. Assume elastic collision (kinetic energy is conserved)

        // 5a. Conservation of linear momentum
        const bool iForcedX = shapes[i]->bForcedInSimFrame[0];
        const bool iForcedY = shapes[i]->bForcedInSimFrame[1];
        const bool iForcedZ = shapes[i]->bForcedInSimFrame[2];
        const bool jForcedX = shapes[j]->bForcedInSimFrame[0];
        const bool jForcedY = shapes[j]->bForcedInSimFrame[1];
        const bool jForcedZ = shapes[j]->bForcedInSimFrame[2];

        // Forced objects are treated as if they had infinite mass
        const Real iInvMassX = iForcedX? 0 : 1/shapes[i]->mass;
        const Real iInvMassY = iForcedY? 0 : 1/shapes[i]->mass;
        const Real iInvMassZ = iForcedZ? 0 : 1/shapes[i]->mass;
        const Real jInvMassX = jForcedX? 0 : 1/shapes[j]->mass;
        const Real jInvMassY = jForcedY? 0 : 1/shapes[j]->mass;
        const Real jInvMassZ = jForcedZ? 0 : 1/shapes[j]->mass;
        const Real meanMassX = 2 / std::max(iInvMassX + jInvMassX, EPS);
        const Real meanMassY = 2 / std::max(iInvMassY + jInvMassY, EPS);
        const Real meanMassZ = 2 / std::max(iInvMassZ + jInvMassZ, EPS);
        const Real FXdt = meanMassX * projVel * NX;
        const Real FYdt = meanMassY * projVel * NY;
        const Real FZdt = meanMassZ * projVel * NZ;

        // Forced obstacles are not updated because InvMass = 0
        shapes[i]->transVel[0] += FXdt * iInvMassX;
        shapes[i]->transVel[1] += FYdt * iInvMassY;
        shapes[i]->transVel[2] += FZdt * iInvMassY;

        // 5a. Conservation of angular momentum
        const Real iCx =shapes[i]->centerOfMass[0];
        const Real iCy =shapes[i]->centerOfMass[1];
        const Real iCz =shapes[i]->centerOfMass[2];

        const bool iForcedA0 = shapes[i]->bBlockRotation[0];
        const bool iForcedA1 = shapes[i]->bBlockRotation[1];
        const bool iForcedA2 = shapes[i]->bBlockRotation[2];

        const Real RcrossF0 = iForcedA0 ? 0 : (CY-iCy) * FZdt - (CZ-iCz) * FYdt;
        const Real RcrossF1 = iForcedA1 ? 0 : (CZ-iCz) * FXdt - (CX-iCx) * FZdt;
        const Real RcrossF2 = iForcedA2 ? 0 : (CX-iCx) * FYdt - (CY-iCy) * FXdt;

        // J* delta angVel = RcrossF, where J is a 3x3 symmetrix matrix.
        // Here we compute J^{-1} explicitly and do the multiplication J^{-1}*RcrossF
        const Real m11 = shapes[i]->J[0];
        const Real m12 = shapes[i]->J[3];
        const Real m13 = shapes[i]->J[4];
        const Real m22 = shapes[i]->J[1];
        const Real m23 = shapes[i]->J[5];
        const Real m33 = shapes[i]->J[2];

        Real a11 = m33*m22 - m23*m23;
        Real a12 = m13*m23 - m33*m12;
        Real a13 = m12*m23 - m13*m22;
        Real a22 = m33*m11 - m13*m13;
        Real a23 = m12*m13 - m11*m23;
        Real a33 = m11*m22 - m12*m12;
        const Real determinant =  1.0/((m11 * a11) + (m12 * a12) + (m13 * a13));
        a11 *= determinant;
        a12 *= determinant;
        a13 *= determinant;
        a22 *= determinant;
        a23 *= determinant;
        a33 *= determinant;

        shapes[i]->angVel[0] += a11*RcrossF0 + a12*RcrossF1 + a13*RcrossF2;
        shapes[i]->angVel[0] += a12*RcrossF0 + a22*RcrossF1 + a23*RcrossF2;
        shapes[i]->angVel[0] += a13*RcrossF0 + a23*RcrossF1 + a33*RcrossF2;
    }
}


Penalization::Penalization(SimulationData & s) : Operator(s) {}

void Penalization::operator()(const double dt)
{
  if(sim.obstacle_vector->nObstacles() == 0) return;

  preventCollidingObstacles();

  sim.startProfiler("Penalization");
  std::vector<cubism::BlockInfo>& vInfo = sim.vInfo();
  #pragma omp parallel
  { // each thread needs to call its own non-const operator() function
    if(sim.bImplicitPenalization)
    {
      KernelPenalization<1> K(dt, sim.lambda, sim.obstacle_vector);
      #pragma omp for schedule(dynamic, 1)
      for (size_t i = 0; i < vInfo.size(); ++i) K(vInfo[i]);
    }
    else
    {
      KernelPenalization<0> K(dt, sim.lambda, sim.obstacle_vector);
      #pragma omp for schedule(dynamic, 1)
      for (size_t i = 0; i < vInfo.size(); ++i) K(vInfo[i]);
    }
  }

  ObstacleVisitor*K = new KernelFinalizePenalizationForce(sim.grid);
  sim.obstacle_vector->Accept(K); // accept you son of a french cow
  delete K;

  sim.stopProfiler();
  check("Penalization");
}

CubismUP_3D_NAMESPACE_END
