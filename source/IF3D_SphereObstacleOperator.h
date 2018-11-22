//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#pragma once

#include "IF3D_ObstacleOperator.h"

class IF3D_SphereObstacleOperator: public IF3D_ObstacleOperator
{
  const double radius;
  //special case: startup with unif accel to umax in tmax, and then decel to 0
  bool accel_decel = false;
  double umax = 0, tmax = 1, tOld = 0;

public:

  IF3D_SphereObstacleOperator(FluidGridMPI *g, ArgumentParser &p, const Real *u);
  void create(const int step_id,const double time, const double dt, const Real *Uinf) override;
  void finalize(const int step_id,const double time, const double dt, const Real *Uinf) override;
  void computeVelocities(const Real* Uinf) override;
};
