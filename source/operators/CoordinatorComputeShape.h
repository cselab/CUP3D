//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#ifndef CubismUP_3D_CoordinatorComputeShape_h
#define CubismUP_3D_CoordinatorComputeShape_h

#include "GenericCoordinator.h"
#include "../obstacles/IF3D_ObstacleVector.h"

struct VelocityObstacleVisitor : public ObstacleVisitor
{
  FluidGridMPI * grid;
  const Real * const uInf;
  int * const nSum;
  double * const uSum;
  std::vector<BlockInfo> vInfo;
  bool finalize = false;

  VelocityObstacleVisitor(FluidGridMPI* _grid, const Real*const _uInf,
    int*const _nSum, double*const _uSum) : grid(_grid), uInf(_uInf), nSum(_nSum), uSum(_uSum)
  {
    vInfo = grid->getBlocksInfo();
  }

  void visit(IF3D_ObstacleOperator* const obstacle)
  {
    if (not finalize)
    {
      const auto &bFixFrameOfRef = obstacle->bFixFrameOfRef;
      const Real dummy[3] = { 0.0, 0.0, 0.0 };
      obstacle->computeVelocities(dummy); // compute velocities with zero uinf
      double povU[3];
      obstacle->getTranslationVelocity(povU);

      if (bFixFrameOfRef[0]) { (nSum[0])++; uSum[0] -= povU[0]; }
      if (bFixFrameOfRef[1]) { (nSum[1])++; uSum[1] -= povU[1]; }
      if (bFixFrameOfRef[2]) { (nSum[2])++; uSum[2] -= povU[2]; }
    }
    else
    {
      double obstU[3];
      obstacle->getTranslationVelocity(obstU);
        obstU[0] += uInf[0]; //set obstacle speed to zero
        obstU[1] += uInf[1];
        obstU[2] += uInf[2];
      obstacle->setTranslationVelocity(obstU);
    }
  }
};

class CoordinatorComputeShape : public GenericCoordinator
{

 public:
  CoordinatorComputeShape(SimulationData & s) : GenericCoordinator(s) {}

  void operator()(const double dt)
  {
    check("shape - start");
    sim.obstacle_vector->update(sim.step, sim.time, dt, sim.uinf);

    //each obstacle does chi = max(chi,obstacleblock->chi)
    #pragma omp parallel
    {
      #pragma omp for schedule(static)
      for (int i = 0; i < (int)vInfo.size(); ++i) {
        BlockInfo info = vInfo[i];
        FluidBlock& b = *(FluidBlock*)info.ptrBlock;

        for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
        for(int iy=0; iy<FluidBlock::sizeY; ++iy)
        for(int ix=0; ix<FluidBlock::sizeX; ++ix) b(ix,iy,iz).chi = 0;
      }
    }

    sim.obstacle_vector->create(sim.step, sim.time, dt, sim.uinf);

    {
     int nSum[3] = {0,0,0};
     double uSum[3] = {0,0,0};
     ObstacleVisitor*VIS = new VelocityObstacleVisitor(grid,sim.uinf,nSum,uSum);
     sim.obstacle_vector->Accept(VIS);//accept you son of a french cow
     if(nSum[0]) sim.uinf[0] = uSum[0]/nSum[0];
     if(nSum[1]) sim.uinf[1] = uSum[1]/nSum[1];
     if(nSum[2]) sim.uinf[2] = uSum[2]/nSum[2];
     //printf("Old Uinf %g %g %g\n",uInf[0],uInf[1],uInf[2]);
     VelocityObstacleVisitor*VVIS = static_cast<VelocityObstacleVisitor*>(VIS);
     VVIS->finalize = true;
     sim.obstacle_vector->Accept(VIS);//accept you son of a french cow
     //if(rank == 0) if(nSum[0] || nSum[1] || nSum[2])
     //  printf("New Uinf %g %g %g (from %d %d %d)\n",
     //  uInf[0],uInf[1],uInf[2],nSum[0],nSum[1],nSum[2]);
     delete VIS;
    }
    check("shape - end");
  }

  std::string getName() {
    return "ComputeShape";
  }
};


class CoordinatorComputeForces : public GenericCoordinator
{
 public:
  CoordinatorComputeForces(SimulationData & s) : GenericCoordinator(s) {}

  void operator()(const double dt)
  {
    check((std::string)"obst. forces - start");
    const Real time = sim.time, nu = sim.nu; const auto step = sim.step;
    sim.obstacle_vector->computeForces(step, time, dt, sim.uinf, nu, sim.bDump);
    check((std::string)"obst. forces - end");
  }

  std::string getName() { return "ComputeForces"; }
};


class CoordinatorComputeDiagnostics : public GenericCoordinator
{
 public:
  CoordinatorComputeDiagnostics(SimulationData & s) : GenericCoordinator(s) {}

  void operator()(const double dt)
  {
    check((std::string)"obst. diagnostics - start");
    const Real time = sim.time, lambda = sim.lambda;
    sim.obstacle_vector->computeDiagnostics(sim.step, time, sim.uinf, lambda);
    check((std::string)"obst. diagnostics - end");
  }

  std::string getName() { return "ComputeObstDiag"; }
};

#endif
