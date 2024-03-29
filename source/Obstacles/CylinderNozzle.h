//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//

#pragma once

#include "Obstacle.h"
#include "Cylinder.h"
#include "extra/Schedulers.h"

CubismUP_3D_NAMESPACE_BEGIN

class CylinderNozzle : public Cylinder
{
  std::vector<Real> actuators_prev_value;
  std::vector<Real> actuators_next_value;
  const int Nactuators;
  const Real actuator_theta;
  Real fx_integral = 0;
  std::vector<double>   action_taken;
  std::vector<double> t_action_taken;
  int curr_idx = 0;
  std::vector < Schedulers::ParameterSchedulerScalar > actuatorSchedulers;
  Real t_change = 0;
  const Real regularizer;
  const Real ccoef;
public:
  CylinderNozzle(SimulationData&s, cubism::ArgumentParser &p);
  void finalize() override;
  void act(std::vector<Real> action, const int agentID);
  Real reward(const int agentID);
  std::vector<Real> state(const int agentID);
  std::vector<Real> actuators;
  void computeVelocities() override
  {
       Obstacle::computeVelocities();
       constexpr Real t1 = 0.25;
       constexpr Real t2 = 0.50;
       angVel[2] = (sim.time > t1 && sim.time < t2) ? transVel[0]*2*radius*sin(2*M_PI*(sim.time-t1)/(t2-t1)) : 0.0;
  }
};

CubismUP_3D_NAMESPACE_END
