//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#ifndef CubismUP_3D_PenalizationObstacleVisitor_h
#define CubismUP_3D_PenalizationObstacleVisitor_h

#include "obstacles/ObstacleVector.h"
#include <cmath>

CubismUP_3D_NAMESPACE_BEGIN
using namespace cubism;

struct PenalizationObstacleVisitor : public ObstacleVisitor
{
  FluidGridMPI * const grid;
  const double dt;
  const Real * const uInf;
  const std::vector<BlockInfo>& vInfo = grid->getBlocksInfo();

  PenalizationObstacleVisitor(FluidGridMPI*g, const double _dt,
    const Real*const u) : grid(g), dt(_dt), uInf(u) { }

  void visit(Obstacle* const obstacle)
  {
    //using CHI_MAT = Real[CUP_BLOCK_SIZE][CUP_BLOCK_SIZE][CUP_BLOCK_SIZE];
    using UDEFMAT = Real[CUP_BLOCK_SIZE][CUP_BLOCK_SIZE][CUP_BLOCK_SIZE][3];
    #pragma omp parallel
    {
      const std::vector<ObstacleBlock*>& obstblocks =
        obstacle->getObstacleBlocks();
      double uBody[3], omegaBody[3], centerOfMass[3];
      obstacle->getCenterOfMass(centerOfMass);
      obstacle->getTranslationVelocity(uBody);
      obstacle->getAngularVelocity(omegaBody);
      const size_t Nblocks = vInfo.size();
      #pragma omp for schedule(dynamic)
      for (size_t i = 0; i < Nblocks; ++i)
      {
        const BlockInfo& info = vInfo[i];
        const auto pos = obstblocks[info.blockID];
        if(pos == nullptr) continue;

        FluidBlock& b = *(FluidBlock*)info.ptrBlock;
        UDEFMAT & __restrict__ UDEF = pos->udef;

        for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
        for(int iy=0; iy<FluidBlock::sizeY; ++iy)
        for(int ix=0; ix<FluidBlock::sizeX; ++ix)
        {
          Real p[3]; info.pos(p, ix, iy, iz);
          p[0]-=centerOfMass[0]; p[1]-=centerOfMass[1]; p[2]-=centerOfMass[2];
          const Real object_UR[3] = {
              (Real) omegaBody[1]*p[2] - (Real) omegaBody[2]*p[1],
              (Real) omegaBody[2]*p[0] - (Real) omegaBody[0]*p[2],
              (Real) omegaBody[0]*p[1] - (Real) omegaBody[1]*p[0]
          };
          const Real U_TOT[3] = {
              (Real)uBody[0] + object_UR[0] + UDEF[iz][iy][ix][0],
              (Real)uBody[1] + object_UR[1] + UDEF[iz][iy][ix][1],
              (Real)uBody[2] + object_UR[2] + UDEF[iz][iy][ix][2]
          };
          // what if multiple obstacles share a block??
          // let's plus equal and wake up during the night to stress about it
          b(ix,iy,iz).tmpU += U_TOT[0];
          b(ix,iy,iz).tmpV += U_TOT[1];
          b(ix,iy,iz).tmpW += U_TOT[2];
        }
      }
    }
  }
};

struct InitialPenalization : public ObstacleVisitor
{
  FluidGridMPI * const grid;
  const double dt;
  const Real * const uInf;
  const std::vector<BlockInfo>& vInfo = grid->getBlocksInfo();

  InitialPenalization(FluidGridMPI*g, const double _dt,
    const Real*const u) : grid(g), dt(_dt), uInf(u) { }

  void visit(Obstacle* const obstacle)
  {
    using CHI_MAT = Real[CUP_BLOCK_SIZE][CUP_BLOCK_SIZE][CUP_BLOCK_SIZE];
    using UDEFMAT = Real[CUP_BLOCK_SIZE][CUP_BLOCK_SIZE][CUP_BLOCK_SIZE][3];
    #pragma omp parallel
    {
      const std::vector<ObstacleBlock*>& obstblocks =
        obstacle->getObstacleBlocks();
      double uBody[3], omegaBody[3], centerOfMass[3];
      obstacle->getCenterOfMass(centerOfMass);
      obstacle->getTranslationVelocity(uBody);
      obstacle->getAngularVelocity(omegaBody);
      const size_t Nblocks = vInfo.size();
      #pragma omp for schedule(dynamic)
      for (size_t i = 0; i < Nblocks; ++i)
      {
        const BlockInfo& info = vInfo[i];
        const auto pos = obstblocks[info.blockID];
        if(pos == nullptr) continue;

        FluidBlock& b = *(FluidBlock*)info.ptrBlock;
        CHI_MAT & __restrict__ CHI = pos->chi;
        UDEFMAT & __restrict__ UDEF = pos->udef;

        for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
        for(int iy=0; iy<FluidBlock::sizeY; ++iy)
        for(int ix=0; ix<FluidBlock::sizeX; ++ix)
        {
          Real p[3]; info.pos(p, ix, iy, iz);
          p[0]-=centerOfMass[0]; p[1]-=centerOfMass[1]; p[2]-=centerOfMass[2];
          const Real object_UR[3] = {
              (Real) omegaBody[1]*p[2] - (Real) omegaBody[2]*p[1],
              (Real) omegaBody[2]*p[0] - (Real) omegaBody[0]*p[2],
              (Real) omegaBody[0]*p[1] - (Real) omegaBody[1]*p[0]
          };
          const Real U_TOT[3] = {
              (Real)uBody[0] + object_UR[0] + UDEF[iz][iy][ix][0],
              (Real)uBody[1] + object_UR[1] + UDEF[iz][iy][ix][1],
              (Real)uBody[2] + object_UR[2] + UDEF[iz][iy][ix][2]
          };
          // what if multiple obstacles share a block??
          // let's plus equal and wake up during the night to stress about it
          b(ix,iy,iz).u += CHI[iz][iy][ix] * ( U_TOT[0] - b(ix,iy,iz).u );
          b(ix,iy,iz).v += CHI[iz][iy][ix] * ( U_TOT[1] - b(ix,iy,iz).v );
          b(ix,iy,iz).w += CHI[iz][iy][ix] * ( U_TOT[2] - b(ix,iy,iz).w );
        }
      }
    }
  }
};

CubismUP_3D_NAMESPACE_END
#endif
