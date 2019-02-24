//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#include "ObstacleManagement.h"
#include "../obstacles/IF3D_ObstacleVector.h"

struct VelocityObstacleVisitor : public ObstacleVisitor
{
  const double dt, lambda;
  int * const nSum;
  double * const uSum;

  VelocityObstacleVisitor(SimulationData& s, int*const nS, double*const uS) :
   dt(s.dt), lambda(s.lambda), nSum(nS), uSum(uS) { }

  void visit(IF3D_ObstacleOperator* const obstacle)
  {
    const auto &bFixFrameOfRef = obstacle->bFixFrameOfRef;
    obstacle->computeVelocities(dt, lambda);
    double povU[3];
    obstacle->getTranslationVelocity(povU);

    if (bFixFrameOfRef[0]) { (nSum[0])++; uSum[0] -= povU[0]; }
    if (bFixFrameOfRef[1]) { (nSum[1])++; uSum[1] -= povU[1]; }
    if (bFixFrameOfRef[2]) { (nSum[2])++; uSum[2] -= povU[2]; }
  }
};

void CreateObstacles::operator()(const double dt)
{
  sim.startProfiler("Obst. Reset");
  #pragma omp parallel for schedule(static)
  for (size_t i = 0; i < vInfo.size(); ++i)
  {
    FluidBlock& b = *(FluidBlock*)vInfo[i].ptrBlock;
    for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
    for(int iy=0; iy<FluidBlock::sizeY; ++iy)
    for(int ix=0; ix<FluidBlock::sizeX; ++ix) b(ix,iy,iz).chi = 0;
  }
  sim.stopProfiler();

  sim.startProfiler("Obst. Create");
  sim.obstacle_vector->create(sim.step, sim.time, dt, sim.uinf);
  sim.stopProfiler();
  check("obst. create");
}

void ComputeForces::operator()(const double dt)
{
  sim.startProfiler("Obst. Forces");
  const Real time = sim.time, nu = sim.nu; const auto step = sim.step;
  sim.obstacle_vector->computeForces(step, time, dt, nu, sim.bDump);
  sim.stopProfiler();
  check("obst. forces");
}

void UpdateObstacles::operator()(const double dt)
{
  int nSum[3] = {0,0,0};
  double uSum[3] = {0,0,0};

  sim.startProfiler("Obst. Update");
  ObstacleVisitor* VIS = new VelocityObstacleVisitor(sim, nSum, uSum);
  sim.obstacle_vector->Accept(VIS); // accept you son of a french cow
  if( nSum[0] > 0 ) sim.uinf[0] = uSum[0] / nSum[0];
  if( nSum[1] > 0 ) sim.uinf[1] = uSum[1] / nSum[1];
  if( nSum[2] > 0 ) sim.uinf[2] = uSum[2] / nSum[2];
  //if(rank == 0) if(nSum[0] || nSum[1] || nSum[2])
  //  printf("New Uinf %g %g %g (from %d %d %d)\n",
  //  uInf[0],uInf[1],uInf[2],nSum[0],nSum[1],nSum[2]);
  delete VIS;

  sim.obstacle_vector->update(sim.step, sim.time, dt, sim.uinf);

  sim.stopProfiler();
  check("obst. update");
}
