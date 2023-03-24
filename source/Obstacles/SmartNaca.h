//
//  CubismUP_3D
//  Copyright (c) 2023 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//

#pragma once

#include "Naca.h"
#include "extra/Schedulers.h"

CubismUP_3D_NAMESPACE_BEGIN

class SmartNaca: public Naca
{
  std::vector<Real> actuators;
  std::vector<Real> actuators_prev_value;
  std::vector<Real> actuators_next_value;
  const int Nactuators;
  const Real actuator_ds;
  Real fx_integral = 0;
  std::vector < Schedulers::ParameterSchedulerScalar > actuatorSchedulers;
  Real t_change = 0;
  const Real thickness;

 public:
  SmartNaca(SimulationData&s, cubism::ArgumentParser&p);
  void finalize() override;
  void act( std::vector<Real> action, const int agentID);
  Real reward(const int agentID);
  std::vector<Real> state(const int agentID);
};

CubismUP_3D_NAMESPACE_END
