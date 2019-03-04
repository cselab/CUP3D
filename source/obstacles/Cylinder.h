//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#ifndef CubismUP_3D_Cylinder_h
#define CubismUP_3D_Cylinder_h

#include "obstacles/Obstacle.h"

CubismUP_3D_NAMESPACE_BEGIN

class Cylinder : public Obstacle
{
  const double radius;
  const double halflength;
  std::string section = "circular"; // or whatever
public:
  Cylinder(SimulationData&s, ArgumentParser &p);
  Cylinder(SimulationData&s, ObstacleArguments &args,
                                 double radius, double halflength);

  void _init(void);
  void create() override;
  void finalize() override;
};

CubismUP_3D_NAMESPACE_END
#endif // CubismUP_3D_Cylinder_h
