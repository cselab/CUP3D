//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#ifndef CubismUP_2D_OperatorComputeForces_h
#define CubismUP_2D_OperatorComputeForces_h

#include "operators/Operator.h"
#include "obstacles/IF3D_ObstacleVector.h"

class ComputeForces : public Operator
{
 public:
  ComputeForces(SimulationData & s) : Operator(s) {}

  void operator()(const double dt);

  std::string getName() { return "ComputeForces"; }
};

#endif
