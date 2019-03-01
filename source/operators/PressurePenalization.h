//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#ifndef CubismUP_3D_PressurePenalization_h
#define CubismUP_3D_PressurePenalization_h

#include "operators/Operator.h"

class PressurePenalization : public Operator
{
 protected:
  PoissonSolver * pressureSolver;

 public:
  PressurePenalization(SimulationData & s);

  void operator()(const double dt);

  std::string getName() { return "PressurePenalization"; }
};
#endif
