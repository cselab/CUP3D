//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#ifndef CubismUP_3D_CoordinatorFadeOut_h
#define CubismUP_3D_CoordinatorFadeOut_h

#include "Operator.h"

class FadeOut : public Operator
{
 public:
  FadeOut(SimulationData & s) : Operator(s) { }
  void operator()(const double dt);
  std::string getName() { return "FadeOut"; }
};

#endif
