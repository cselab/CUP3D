//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#ifndef CubismUP_3D_ObstaclesUpdate_h
#define CubismUP_3D_ObstaclesUpdate_h

#include "Operator.h"

CubismUP_3D_NAMESPACE_BEGIN

class UpdateObstacles : public Operator
{
 public:
  UpdateObstacles(SimulationData & s) : Operator(s) {}

  void operator()(const Real dt);

  std::string getName() { return "UpdateObstacles Vel"; }
};

CubismUP_3D_NAMESPACE_END
#endif // CubismUP_3D_ObstaclesUpdate_h
