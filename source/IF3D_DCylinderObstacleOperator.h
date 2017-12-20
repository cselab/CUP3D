//
//  CubismUP_3D
//
//  Written by Guido Novati ( novatig@ethz.ch ).
//  This file started as an extension of code written by Wim van Rees
//  Copyright (c) 2017 ETHZ. All rights reserved.
//

#pragma once

#include "IF3D_ObstacleOperator.h"

#include <cmath>

class IF3D_DCylinderObstacleOperator: public IF3D_ObstacleOperator
{
  double radius, halflength;

public:

 IF3D_DCylinderObstacleOperator(FluidGridMPI*g, ArgumentParser&p, const Real*const u) //const Real radius, const double position[3], const Real smoothing_length=-1):
  : IF3D_ObstacleOperator(g, p, u)//, radius(radius), smoothing_length(smoothing_length)
  {
      _parseArguments(p);
  }

  void create(const int step_id,const double time, const double dt, const Real *Uinf) override;

  void _parseArguments(ArgumentParser & parser) override;
};
