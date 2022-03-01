//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#include "ObstacleVector.h"

#include <sstream>

CubismUP_3D_NAMESPACE_BEGIN
using namespace cubism;

std::vector<std::array<int, 2>> ObstacleVector::collidingObstacles()
{
  std::set<std::array<int, 2>> colliding; //IDs of colliding obstacles
  //vector containing pointers to defBLock maps:
  std::vector<std::vector<ObstacleBlock*>*> obstBlocks(obstacles.size());
  int ID = 0;
  std::vector<int> IDs;
  for(const auto & obstacle_ptr : obstacles) {
      IDs.push_back(ID);
      obstBlocks[ID] = obstacle_ptr->getObstacleBlocksPtr();
      assert(obstBlocks[ID] != nullptr);
      ID++;
  }

  for(int i=1; i<ID; i++)
  for(int j=0; j<i; j++) {
    const auto& y = * obstBlocks[j];
    const auto& x = * obstBlocks[i];
    assert(x.size() == y.size());
    for (size_t k=0; k<x.size(); k++) {
      if(x[k] not_eq nullptr && y[k] not_eq nullptr) {
        std::array<int,2> hit = {IDs[i],IDs[j]};
        colliding.insert(hit); //it's a set: only unique pairs are inserted
      }
    }
  }
  return std::vector<std::array<int,2>>(colliding.begin(), colliding.end());
}

void ObstacleVector::update()
{
    for(const auto & obstacle_ptr : obstacles)
        obstacle_ptr->update();
}

void ObstacleVector::create()
{
  for(const auto & obstacle_ptr : obstacles)
    obstacle_ptr->create();
}

void ObstacleVector::finalize()
{
  for(const auto & obstacle_ptr : obstacles)
    obstacle_ptr->finalize();
}

void ObstacleVector::computeVelocities()
{
  for(const auto & obstacle_ptr : obstacles)
    obstacle_ptr->computeVelocities();
}

void ObstacleVector::computeForces()
{
  for(const auto & obstacle_ptr : obstacles)
    obstacle_ptr->computeForces();
}

void ObstacleVector::save(std::string filename)
{
  int cntr = 0;
  for(const auto & obstacle_ptr : obstacles) {
    std::stringstream ssR;
    ssR<<filename<<"_"<<cntr;
      obstacle_ptr->save(ssR.str());
      cntr++;
  }
}

void ObstacleVector::restart(std::string filename)
{
    int cntr = 0;
    for(const auto & obstacle_ptr : obstacles) {
      std::stringstream ssR;
      ssR<<filename<<"_"<<cntr;
        obstacle_ptr->restart(ssR.str());
        cntr++;
    }
}

void ObstacleVector::Accept(ObstacleVisitor * visitor)
{
  for(size_t i=0;i<obstacles.size();++i)
    obstacles[i]->Accept(visitor);
}

CubismUP_3D_NAMESPACE_END
