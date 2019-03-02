//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#pragma once

#include "obstacles/IF3D_ObstacleOperator.h"

class IF3D_EllipsoidObstacleOperator: public IF3D_ObstacleOperator
{
  const double radius;
  Real e0, e1, e2;
  //special case: startup with unif accel to umax in tmax, and then decel to 0
  bool accel_decel = false;
  double umax = 0, tmax = 1;

public:

  IF3D_EllipsoidObstacleOperator(SimulationData&s,ArgumentParser&p);

  void create() override;
  void finalize() override;
  void computeVelocities() override;
};
