//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#pragma once

#include "Operator.h"

CubismUP_3D_NAMESPACE_BEGIN

class AdvectionDiffusion : public Operator
{
  std::vector<Real> vOld;
public:
  AdvectionDiffusion(SimulationData & s) : Operator(s) { }

  ~AdvectionDiffusion() { }

  void operator()(const Real dt);

  std::string getName() { return "AdvectionDiffusion"; }
};

CubismUP_3D_NAMESPACE_END
