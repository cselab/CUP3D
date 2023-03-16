//
//  Cubism3D
//  Copyright (c) 2023 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//

#pragma once

#include "../SimulationData.h"
#include "Operator.h"

CubismUP_3D_NAMESPACE_BEGIN

class ExternalForcing : public Operator
{
 public:
  ExternalForcing(SimulationData & s) : Operator(s) {}

  void operator()(const double dt);

  std::string getName() { return "ExternalForcing"; }
};

CubismUP_3D_NAMESPACE_END
