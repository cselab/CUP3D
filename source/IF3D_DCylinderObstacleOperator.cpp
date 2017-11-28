//
//  CubismUP_3D
//
//  Written by Guido Novati ( novatig@ethz.ch ).
//  This file started as an extension of code written by Wim van Rees
//  Copyright (c) 2017 ETHZ. All rights reserved.
//

#include <limits>
#include <cstring>

#include "IF3D_DCylinderObstacleOperator.h"
#include "IF3D_ObstacleLibrary.h"

void IF3D_DCylinderObstacleOperator::create(const int step_id,const double time, const double dt, const Real *Uinf)
{
  for(auto & entry : obstacleBlocks) delete entry.second;
  obstacleBlocks.clear();

  DCylinderObstacle::FillBlocks kernel(radius, halflength, vInfo[0].h_gridpoint, position);
  for(int i=0; i<vInfo.size(); i++) {
    BlockInfo info = vInfo[i];
    //const auto pos = obstacleBlocks.find(info.blockID);
    if(kernel._is_touching(info)) { //position of sphere + radius + 2*h safety
      assert(obstacleBlocks.find(info.blockID) == obstacleBlocks.end());
      obstacleBlocks[info.blockID] = new ObstacleBlock;
      obstacleBlocks[info.blockID]->clear(); //memset 0
    }
  }

  //this writes the chi field, therefore moved to finalize
  #pragma omp parallel
  {
    DCylinderObstacle::FillBlocks fill(radius, halflength, vInfo[0].h_gridpoint, position);
    const int tid = omp_get_thread_num();

    #pragma omp for schedule(dynamic)
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

void IF3D_DCylinderObstacleOperator::setTranslationVelocity(double UT[3])
{}

void IF3D_DCylinderObstacleOperator::_parseArguments(ArgumentParser & parser)
{
  //obstacleop parses x,y,z,quats and length!
  IF3D_ObstacleOperator::_parseArguments(parser);
  parser.set_strict_mode();
  parser.unset_strict_mode();
  const double xvel = parser("-xvel").asDouble(0.);
  const double yvel = parser("-yvel").asDouble(0.);
  const double zvel = parser("-zvel").asDouble(0.);

  this->transVel[0] = xvel;
  this->transVel[1] = yvel;
  this->transVel[2] = zvel;
  halflength = ext_Z/2;
  radius = .5*length;
  printf("Created IF3D_DCylinderObstacleOperator with radius %f\n", radius);
}