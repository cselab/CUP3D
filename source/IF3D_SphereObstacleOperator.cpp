//
//  CubismUP_3D
//
//  Written by Guido Novati ( novatig@ethz.ch ).
//  This file started as an extension of code written by Wim van Rees
//  Copyright (c) 2017 ETHZ. All rights reserved.
//

#include <limits>
#include <cstring>

#include "IF3D_SphereObstacleOperator.h"
#include "IF3D_ObstacleLibrary.h"

void IF3D_SphereObstacleOperator::create(const int step_id, const double time, const double dt, const Real *Uinf)
{
  for(auto & entry : obstacleBlocks) delete entry.second;
  obstacleBlocks.clear();

  SphereObstacle::FillBlocks kernel(radius,vInfo[0].h_gridpoint,position);
  for(int i=0; i<vInfo.size(); i++) {
    BlockInfo info = vInfo[i];
    //const auto pos = obstacleBlocks.find(info.blockID);
    if(kernel._is_touching(info)) { //position of sphere + radius + 2*h safety
      assert(obstacleBlocks.find(info.blockID) == obstacleBlocks.end());
      obstacleBlocks[info.blockID] = new ObstacleBlock();
      obstacleBlocks[info.blockID]->clear(); //memset 0
    }
  }
  tOld = time;
}

void IF3D_SphereObstacleOperator::finalize(const int step_id,const double time, const double dt, const Real *Uinf)
{
  //this writes the chi field, therefore moved to finalize
  #pragma omp parallel
  {
    SphereObstacle::FillBlocks fill(radius, vInfo[0].h_gridpoint, position);
    const int tid = omp_get_thread_num();

    #pragma omp for schedule(static)
    for(int i=0; i<vInfo.size(); i++) {
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
    if(tOld<tmax) transVel[0] = umax*tOld/tmax;
    else if (tOld<2*tmax) transVel[0] = umax*(2*tmax-tOld)/tmax;
    else transVel[0] = 0;
  }
}
