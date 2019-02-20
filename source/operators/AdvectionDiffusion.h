//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#ifndef CubismUP_3D_AdvectionDiffusion_h
#define CubismUP_3D_AdvectionDiffusion_h

#include "Operator.h"

class AdvectionDiffusion : public Operator
{
public:
  AdvectionDiffusion(SimulationData & s) : Operator(s) { }

  ~AdvectionDiffusion() { }

  void operator()(const double dt);

  std::string getName() { return "AdvectionDiffusion"; }
};

#endif
