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

  IF3D_SphereObstacleOperator(FluidGridMPI*g, ArgumentParser&p, const Real*const u)
  : IF3D_ObstacleOperator(g, p, u), radius(0.5*length)
  {
    accel_decel = p("-accel").asBool(false);
    if(accel_decel) {
      if(not bForcedInSimFrame[0]) {
        printf("Warning: sphere was not set to be forced in x-dir, yet the accel_decel pattern is active.\n");
      }
      umax = p("-xvel").asDouble(0.0);
      tmax = p("-T").asDouble(1.);
    }
  }

  void create(const int step_id,const double time, const double dt, const Real *Uinf) override;
  void finalize(const int step_id,const double time, const double dt, const Real *Uinf) override;
  void computeVelocities(const Real* Uinf) override;
};
