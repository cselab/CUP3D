//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#ifndef CubismUP_3D_PressureRHS_h
#define CubismUP_3D_PressureRHS_h

#include "Operator.h"

class PressureRHS : public Operator
{
 public:
  PressureRHS(SimulationData & s) : Operator(s) {}

  void operator()(const double dt);

  std::string getName() { return "PressureRHS"; }
};
#endif
