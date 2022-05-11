//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#ifndef CubismUP_3D_Cylinder_h
#define CubismUP_3D_Cylinder_h

#include "Obstacle.h"

CubismUP_3D_NAMESPACE_BEGIN

class Cylinder : public Obstacle
{
  const Real radius;
  const Real halflength;
  std::string section = "circular"; // or whatever
  Real umax = 0;
  Real vmax = 0;
  Real wmax = 0;
  Real tmax = 1;
  //special case: startup with unif accel to umax in tmax, and then decel to 0
  bool accel = false;

public:
  Cylinder(SimulationData&s, cubism::ArgumentParser &p);
  Cylinder(SimulationData&s, ObstacleArguments &args,
                                 Real radius, Real halflength);

  void _init(void);
  void create() override;
  void finalize() override;
  void computeVelocities() override;
};

CubismUP_3D_NAMESPACE_END
#endif // CubismUP_3D_Cylinder_h
