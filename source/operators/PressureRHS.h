//
//  CubismUP_3D
//  Copyright (c) 2020 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Michalis Chatzimanolakis (michaich@ethz.ch).
//

#pragma once

#include "Operator.h"

#include "Cubism/FluxCorrectionMPI.h"

CubismUP_3D_NAMESPACE_BEGIN
#define PENAL_THEN_PRES

class PressureRHS : public Operator
{
  PenalizationGridMPI * penalizationGrid = nullptr;
 public:
  PressureRHS(SimulationData & s);
  ~PressureRHS();

  void operator()(const double dt);

  std::string getName() { return "PressureRHS"; }
};

CubismUP_3D_NAMESPACE_END