//
//  CubismUP_3D
//  Copyright (c) 2020 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Michalis Chatzimanolakis (michaich@ethz.ch).
//

#pragma once

#include <memory>
#include "Operator.h"
#include "../poisson/PoissonSolverBase.h"
#include "../Obstacles/ObstacleVector.h"
#include "../poisson/PoissonSolverAMR.h"

CubismUP_3D_NAMESPACE_BEGIN

class PoissonSolverBase;

class PressureProjection : public Operator
{
 protected:
  // Alias of sim.pressureSolver.
 std::shared_ptr<PoissonSolverBase> pressureSolver;
 std::vector<Real> pOld;

 public:
  PressureProjection(SimulationData & s); 
  ~PressureProjection() = default;

  void operator()(Real dt) override;

  std::string getName() { return "PressureProjection"; }
};

CubismUP_3D_NAMESPACE_END
