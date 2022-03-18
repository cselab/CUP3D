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

class PoissonSolverAMR;

class PressureProjection : public Operator
{
 protected:
  // Alias of sim.pressureSolver.
  PoissonSolverAMR * pressureSolver;

 public:
  PressureProjection(SimulationData & s); 
  ~PressureProjection();

  void operator()(Real dt) override;

  std::string getName() { return "PressureProjection"; }
};

CubismUP_3D_NAMESPACE_END
