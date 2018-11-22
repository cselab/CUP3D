//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#include "IF3D_DCylinderObstacleOperator.h"

#include "Cubism/ArgumentParser.h"
#include "IF3D_ObstacleLibrary.h"

IF3D_DCylinderObstacleOperator::IF3D_DCylinderObstacleOperator(
    FluidGridMPI * const g,
    ArgumentParser &p,
    const Real * const u)
    : IF3D_ObstacleOperator(g, p, u),
      radius(.5 * length),
      halflength(p("-halflength").asDouble(.5 * ext_Z))
{
  printf("Created IF3D_DCylinderObstacleOperator with radius %f and halflength %f\n", radius, halflength);

  // D-cyl can float around the domain, but does not support rotation. TODO
  bBlockRotation[0] = true;
  bBlockRotation[1] = true;
  bBlockRotation[2] = true;
}


void IF3D_DCylinderObstacleOperator::create(const int step_id,const double time, const double dt, const Real *Uinf)
{
  for(auto & entry : obstacleBlocks) delete entry.second;
  obstacleBlocks.clear();

  DCylinderObstacle::FillBlocks kernel(radius, halflength, vInfo[0].h_gridpoint, position);
  for(size_t i=0; i<vInfo.size(); i++) {
    BlockInfo info = vInfo[i];
    //const auto pos = obstacleBlocks.find(info.blockID);
    if(kernel._is_touching(info)) { //position of sphere + radius + 2*h safety
      assert(obstacleBlocks.find(info.blockID) == obstacleBlocks.end());
      obstacleBlocks[info.blockID] = new ObstacleBlock();
      obstacleBlocks[info.blockID]->clear(); //memset 0
    }
  }

  //this writes the chi field, therefore moved to finalize
  #pragma omp parallel
  {
    DCylinderObstacle::FillBlocks fill(radius, halflength, vInfo[0].h_gridpoint, position);

    #pragma omp for schedule(dynamic)
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
