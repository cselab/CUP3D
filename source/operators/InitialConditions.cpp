//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#include "operators/InitialConditions.h"
#include "obstacles/ObstacleVector.h"

CubismUP_3D_NAMESPACE_BEGIN
using namespace cubism;

class KernelIC
{
 public:
  KernelIC(const Real u) {}
  ~KernelIC() {}
  void operator()(const BlockInfo& info, FluidBlock& block) const
  {
    for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
    for(int iy=0; iy<FluidBlock::sizeY; ++iy)
    for(int ix=0; ix<FluidBlock::sizeX; ++ix) block(ix,iy,iz).clear();
  }
};

class KernelIC_RT
{
 public:
  KernelIC_RT(const Real rhoS) {}
  ~KernelIC_RT() {}

  void operator()(const BlockInfo& info, FluidBlock& block) const
  {
    const Real d = .25, a = .05;//.25
    for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
    for(int iy=0; iy<FluidBlock::sizeY; ++iy)
    for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
      block(ix,iy,iz).clear();
      Real p[3];
      info.pos(p, ix, iy, iz);
      Real x = p[0] - d*.5;
      Real y = p[1] - d*2.5;
      Real z = p[2] - d*.5;
      //Real eta = -.1*a*cos(2*M_PI*x/d);
      Real eta = -a*.25*std::cos(2*M_PI*x/d)*std::cos(2*M_PI*z/d);
      block(ix,iy,iz).chi = .5+.5*std::tanh((y-eta)/info.h_gridpoint);
    }
  }
};

class KernelIC_taylorGreen
{
  const Real ext[3], uMax;
  const Real a = 2*M_PI / ext[0], b = 2*M_PI / ext[1], c = 2*M_PI / ext[2];
  const Real A = uMax, B = - uMax * ext[1] / ext[0];
 public:
  ~KernelIC_taylorGreen() {}
  KernelIC_taylorGreen(const Real extent[3], const Real U): ext{extent[0],extent[1],extent[2]}, uMax(U) {}
  void operator()(const BlockInfo& info, FluidBlock& block) const
  {
    for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
    for(int iy=0; iy<FluidBlock::sizeY; ++iy)
    for(int ix=0; ix<FluidBlock::sizeX; ++ix)
    {
      block(ix,iy,iz).clear();
      Real p[3]; info.pos(p, ix, iy, iz);
      block(ix,iy,iz).u = A*std::cos(a*p[0])*std::sin(b*p[1])*std::sin(c*p[2]);
      block(ix,iy,iz).v = B*std::sin(a*p[0])*std::cos(b*p[1])*std::sin(c*p[2]);
    }
  }
};

class KernelIC_channel
{
  const int dir;
  const Real ext[3], uMax, H = ext[dir], FAC = 4*uMax/H/H; // FAC = 0.5*G/mu
  //umax =  0.5*G/mu * 0.25*H*H
 public:
  KernelIC_channel(const Real extent[3], const Real U, const int _dir):
    dir(_dir), ext{extent[0],extent[1],extent[2]}, uMax(U) {}
  ~KernelIC_channel() {}

  void operator()(const BlockInfo& info, FluidBlock& block) const
  {
    for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
    for(int iy=0; iy<FluidBlock::sizeY; ++iy)
    for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
      block(ix,iy,iz).clear();
      Real p[3]; info.pos(p, ix, iy, iz);
      block(ix,iy,iz).u = FAC * p[dir] * (H-p[dir]);
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
      const auto& obstblocks = obstacle->getObstacleBlocks();
      const std::array<double,3> centerOfMass = obstacle->getCenterOfMass();
      const std::array<double,3> uBody = obstacle->getTranslationVelocity();
      const std::array<double,3> omegaBody = obstacle->getAngularVelocity();

      #pragma omp for schedule(dynamic)
      for (size_t i = 0; i < vInfo.size(); ++i)
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

void InitialConditions::operator()(const double dt)
{
  if(sim.initCond == "zero") {
    if(sim.verbose) printf("Zero-values initial conditions.\n");
    run(KernelIC(0));
  }
  if(sim.initCond == "taylorGreen") {
    if(sim.verbose) printf("Taylor Green vortex initial conditions.\n");
    run(KernelIC_taylorGreen(sim.extent, sim.uMax_forced));
  }
  if(sim.initCond == "channel")
  {
    if(sim.verbose) printf("Channel flow initial conditions.\n");
    if( sim.BCx_flag == wall ) {
      printf("ERROR: channel flow must be periodic or dirichlet in x.\n");
      abort();
    }
    const bool channelY = sim.BCy_flag==wall, channelZ = sim.BCz_flag==wall;
    if( (channelY && channelZ) or (!channelY && !channelZ) ) {
      printf("ERROR: wrong channel flow BC in y or z.\n");
      abort();
    }
    const int dir = channelY? 1 : 2;
    run(KernelIC_channel(sim.extent, sim.uMax_forced, dir));
  }
  {
    //zero fields, going to contain Udef:
    #pragma omp parallel for schedule(static)
    for(unsigned i=0; i<vInfo.size(); i++)
    {
      FluidBlock& b = *(FluidBlock*)vInfo[i].ptrBlock;
      for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
      for(int iy=0; iy<FluidBlock::sizeY; ++iy)
      for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
        b(ix,iy,iz).tmpU = 0; b(ix,iy,iz).tmpV = 0; b(ix,iy,iz).tmpW = 0;
      }
    }
    //store deformation velocities onto tmp fields:
    ObstacleVisitor* visitor = new InitialPenalization(grid, dt, sim.uinf);
    sim.obstacle_vector->Accept(visitor);
    delete visitor;
  }

  check("InitialConditions");
}

CubismUP_3D_NAMESPACE_END
