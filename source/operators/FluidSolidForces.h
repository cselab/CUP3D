//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//

#pragma once

#include "Operator.h"
#include "../Obstacles/ObstacleVector.h"

CubismUP_3D_NAMESPACE_BEGIN

class ComputeForces : public Operator
{
 public:
  ComputeForces(SimulationData & s) : Operator(s) {}

  void operator()(const Real dt);

  std::string getName() { return "ComputeForces"; }
};

CubismUP_3D_NAMESPACE_END
