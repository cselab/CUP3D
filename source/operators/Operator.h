//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//

#pragma once
#include "../SimulationData.h"

CubismUP_3D_NAMESPACE_BEGIN

class Operator
{
  public:
  SimulationData& sim;
  Operator(SimulationData& s) noexcept : sim(s) {  }
  virtual ~Operator() = default;
  virtual void operator()(Real dt) = 0;
  virtual std::string getName() = 0;
};
CubismUP_3D_NAMESPACE_END
