//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#include "IF3D_SphereObstacleOperator.h"

#include "IF3D_ObstacleLibrary.h"

#include "Cubism/ArgumentParser.h"

IF3D_SphereObstacleOperator::IF3D_SphereObstacleOperator(
    SimulationData& s, ArgumentParser& p )
    : IF3D_ObstacleOperator(s, p), radius(0.5*length)
{
  accel_decel = p("-accel").asBool(false);
  if(accel_decel) {
    if(not bForcedInSimFrame[0]) {
      printf("Warning: sphere was not set to be forced in x-dir, yet the accel_decel pattern is active.\n");
    }
    umax = p("-xvel").asDouble(0.0);
    tmax = p("-T").asDouble(1.);
  }
}


IF3D_SphereObstacleOperator::IF3D_SphereObstacleOperator(
    SimulationData& s,
    ObstacleArguments &args,
    const double R)
    : IF3D_ObstacleOperator(s, args), radius(R) { }


IF3D_SphereObstacleOperator::IF3D_SphereObstacleOperator(
    SimulationData& s,
    ObstacleArguments &args,
    const double R,
    const double _umax,
    const double _tmax)
    : IF3D_ObstacleOperator(s, args), radius(R), umax(_umax), tmax(_tmax) { }


void IF3D_SphereObstacleOperator::create(const int step_id, const double time, const double dt, const Real *Uinf)
{
  for(auto & entry : obstacleBlocks) delete entry.second;
  obstacleBlocks.clear();

  SphereObstacle::FillBlocks kernel(radius,vInfo[0].h_gridpoint,position);
  for(size_t i=0; i<vInfo.size(); i++) {
    BlockInfo info = vInfo[i];
    //const auto pos = obstacleBlocks.find(info.blockID);
    if(kernel._is_touching(info)) { //position of sphere + radius + 2*h safety
      assert(obstacleBlocks.find(info.blockID) == obstacleBlocks.end());
      obstacleBlocks[info.blockID] = new ObstacleBlock();
      obstacleBlocks[info.blockID]->clear(); //memset 0
    }
  }
  _tOld = time;
}

void IF3D_SphereObstacleOperator::finalize(const int step_id,const double time, const double dt, const Real *Uinf)
{
  //this writes the chi field, therefore moved to finalize
  #pragma omp parallel
  {
    SphereObstacle::FillBlocks fill(radius, vInfo[0].h_gridpoint, position);

    #pragma omp for schedule(static)
    for(size_t i=0; i<vInfo.size(); i++) {
      BlockInfo info = vInfo[i];
      auto pos = obstacleBlocks.find(info.blockID);
      if(pos == obstacleBlocks.end()) continue;
      FluidBlock& b = *(FluidBlock*)info.ptrBlock;
      fill(info, b, pos->second);
    }
  }
  for(auto & o : obstacleBlocks) o.second->allocate_surface();
}


void IF3D_SphereObstacleOperator::computeVelocities(const Real* Uinf)
{
  IF3D_ObstacleOperator::computeVelocities(Uinf);

  if(accel_decel) {
    if(_tOld<tmax) transVel[0] = umax*_tOld/tmax;
    else if (_tOld<2*tmax) transVel[0] = umax*(2*tmax-_tOld)/tmax;
    else transVel[0] = 0;
  }
}
