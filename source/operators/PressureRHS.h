//
//  CubismUP_3D
//  Copyright (c) 2020 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Michalis Chatzimanolakis (michaich@ethz.ch).
//

#pragma once

#include "Operator.h"

CubismUP_3D_NAMESPACE_BEGIN

class PressureRHS : public Operator
{
 public:
  PressureRHS(SimulationData & s);
  ~PressureRHS();

  void operator()(Real dt) override;

  std::string getName() { return "PressureRHS"; }
};

CubismUP_3D_NAMESPACE_END
