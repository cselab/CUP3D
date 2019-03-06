//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#ifndef CubismUP_3D_Sphere_h
#define CubismUP_3D_Sphere_h

#include "obstacles/Obstacle.h"

CubismUP_3D_NAMESPACE_BEGIN

class Sphere: public Obstacle
{
  const double radius;
  //special case: startup with unif accel to umax in tmax, and then decel to 0
  bool accel_decel = false, bHemi = false;
  double umax = 0, tmax = 1;

public:

  Sphere(SimulationData&s,cubism::ArgumentParser&p);
  Sphere(SimulationData&s,ObstacleArguments&args,double R);
  Sphere(SimulationData&s,ObstacleArguments&args,double R, double umax, double tmax);

  void create() override;
  void finalize() override;
  void computeVelocities() override;
};

CubismUP_3D_NAMESPACE_END
#endif // CubismUP_3D_Sphere_h
