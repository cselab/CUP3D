//
//  Cubism3D
//  Copyright (c) 2023 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//

#pragma once

#include "Operator.h"

CubismUP_3D_NAMESPACE_BEGIN

class FixMassFlux : public Operator
{
public:
  FixMassFlux(SimulationData &s);

  void operator()(const double dt);

  std::string getName(){ return "FixMassFlux"; }
};

CubismUP_3D_NAMESPACE_END
