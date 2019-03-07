//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch) and Wim van Rees.
//

#ifndef CubismUP_3D_ObstacleFactory_h
#define CubismUP_3D_ObstacleFactory_h

#include "obstacles/Obstacle.h"
#include <memory>

namespace cubism { struct ArgumentParser; }

CubismUP_3D_NAMESPACE_BEGIN

class ObstacleFactory
{
  SimulationData & sim;
  int _getlines(std::string filename);

public:
  typedef std::vector<std::shared_ptr<Obstacle>> VectorType;

  ObstacleFactory(SimulationData & s) : sim(s) { }
  ~ObstacleFactory() {}
  VectorType create(cubism::ArgumentParser &parser);
};

CubismUP_3D_NAMESPACE_END
#endif // CubismUP_3D_ObstacleFactory_h
