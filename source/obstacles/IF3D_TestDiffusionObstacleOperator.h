//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by SV (vermas@fau.edu).
//

#pragma once

#include "IF3D_ObstacleOperator.h"

class IF3D_TestDiffusionObstacleOperator : public IF3D_ObstacleOperator
{
  const double sigma;

public:
  IF3D_TestDiffusionObstacleOperator(SimulationData& s, ArgumentParser &p);
  IF3D_TestDiffusionObstacleOperator(SimulationData& s, ObstacleArguments &args, double sigma);

  void create(const int step_id,const double time, const double dt, const Real *Uinf) override;
};
