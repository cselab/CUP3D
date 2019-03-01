//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#pragma once

#include "obstacles/IF3D_ObstacleOperator.h"

class IF3D_SphereObstacleOperator: public IF3D_ObstacleOperator
{
  const double radius;
  //special case: startup with unif accel to umax in tmax, and then decel to 0
  bool accel_decel = false;
  double umax = 0, tmax = 1;
  double _tOld = 0;

public:

  IF3D_SphereObstacleOperator(SimulationData&s,ArgumentParser&p);
  IF3D_SphereObstacleOperator(SimulationData&s,ObstacleArguments&args,double R);
  IF3D_SphereObstacleOperator(SimulationData&s,ObstacleArguments&args,double R, double umax, double tmax);

  void create(const int step_id,const double time, const double dt, const Real *Uinf) override;
  void finalize(const int step_id,const double time, const double dt, const Real *Uinf) override;
  void computeVelocities(const double dt, const Real lambda) override;
};
