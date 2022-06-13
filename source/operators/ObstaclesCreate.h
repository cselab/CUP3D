//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//

#pragma once

#include "Operator.h"
#include "../Obstacles/ObstacleVector.h"
#include "../Utils/MatArrayMath.h"

CubismUP_3D_NAMESPACE_BEGIN

class CreateObstacles : public Operator
{
 public:
  CreateObstacles(SimulationData & s) : Operator(s) {}

  void operator()(const Real dt);

  std::string getName() { return "CreateObstacles"; }
};

CubismUP_3D_NAMESPACE_END
