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
  const double radius, halflength;

public:

 IF3D_DCylinderObstacleOperator(FluidGridMPI*g, ArgumentParser&p, const Real*const u) //const Real radius, const double position[3], const Real smoothing_length=-1):
  : IF3D_ObstacleOperator(g, p, u), radius(.5*length), halflength(ext_Z/2)
  {
    printf("Created IF3D_DCylinderObstacleOperator with radius %f\n", radius);
    // D-cyl can float around the domain, but does not support rotation. TODO
    bBlockRotation[0] = true;
    bBlockRotation[1] = true;
    bBlockRotation[2] = true;
  }

  void create(const int step_id,const double time, const double dt, const Real *Uinf) override;
};
