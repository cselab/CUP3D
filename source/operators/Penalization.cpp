//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#include "operators/Penalization.h"
#include "obstacles/ObstacleVector.h"

CubismUP_3D_NAMESPACE_BEGIN
using namespace cubism;

using CHIMAT = Real[CUP_BLOCK_SIZE][CUP_BLOCK_SIZE][CUP_BLOCK_SIZE];
using UDEFMAT = Real[CUP_BLOCK_SIZE][CUP_BLOCK_SIZE][CUP_BLOCK_SIZE][3];

struct KernelPenalization : public ObstacleVisitor
{
  const double lamdt;
  ObstacleVector * const obstacle_vector;
  const cubism::BlockInfo * info_ptr = nullptr;

  KernelPenalization(double lambdadt, ObstacleVector* ov) :
    lamdt(lambdadt), obstacle_vector(ov) {}

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

    FluidBlock& b = *(FluidBlock*)info.ptrBlock;
    const std::array<double,3> CM = obstacle->getCenterOfMass();
    const std::array<double,3> vel = obstacle->getTranslationVelocity();
    const std::array<double,3> omega = obstacle->getAngularVelocity();

    CHIMAT & __restrict__ CHI = o->chi;
    UDEFMAT & __restrict__ UDEF = o->udef;

    for(int iz=0; iz<FluidBlock::sizeZ; iz++)
    for(int iy=0; iy<FluidBlock::sizeY; iy++)
    for(int ix=0; ix<FluidBlock::sizeX; ix++)
    {
      // What if multiple obstacles share a block? Do not write udef onto
      // grid if CHI stored on the grid is greater than obst's CHI.
      if(b(ix,iy,iz).chi > CHI[iz][iy][ix]) continue;
      if(CHI[iz][iy][ix] <= 0) continue; // no need to do anything
      double p[3]; info.pos(p, ix, iy, iz);
      p[0] -= CM[0]; p[1] -= CM[1]; p[2] -= CM[2];

      const double X = CHI[iz][iy][ix], U_TOT[3] = {
          vel[0] + omega[1]*p[2] - omega[2]*p[1] + UDEF[iz][iy][ix][0],
          vel[1] + omega[2]*p[0] - omega[0]*p[2] + UDEF[iz][iy][ix][1],
          vel[2] + omega[0]*p[1] - omega[1]*p[0] + UDEF[iz][iy][ix][2]
      };
      // What if two obstacles overlap? Let's plus equal. After all here
      // we are computing divUs, maybe one obstacle has divUs 0. We will
      // need a repulsion term of the velocity at some point in the code.
      b(ix,iy,iz).u = (b(ix,iy,iz).u + X*lamdt * U_TOT[0]) / (1 + X*lamdt);
      b(ix,iy,iz).v = (b(ix,iy,iz).v + X*lamdt * U_TOT[1]) / (1 + X*lamdt);
      b(ix,iy,iz).w = (b(ix,iy,iz).w + X*lamdt * U_TOT[2]) / (1 + X*lamdt);
    }
  }
};

Penalization::Penalization(SimulationData & s) : Operator(s) {}

void Penalization::operator()(const double dt)
{
  sim.startProfiler("Penalization");
  #pragma omp parallel
  { // each thread needs to call its own non-const operator() function
    auto K = KernelPenalization(dt*sim.lambda, sim.obstacle_vector);
    #pragma omp for schedule(dynamic, 1)
    for (size_t i = 0; i < vInfo.size(); ++i) K(vInfo[i]);
  }
  sim.stopProfiler();

  check("Penalization");
}

CubismUP_3D_NAMESPACE_END
