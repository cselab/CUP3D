//
//  CubismUP_3D
//
//  Written by Guido Novati ( novatig@ethz.ch ).
//  This file started as an extension of code written by Wim van Rees
//  Copyright (c) 2017 ETHZ. All rights reserved.
//

#ifndef __IncompressibleFluids3D__IF3D_ForcedSphereObstacleOperator__
#define __IncompressibleFluids3D__IF3D_ForcedSphereObstacleOperator__

#include "IF3D_SphereObstacleOperator.h"

#include <cmath>

class IF3D_ForcedSphereObstacleOperator: public IF3D_SphereObstacleOperator
{
  bool accel_decel = false;
  double umax = 0, tmax = 1;
public:

  IF3D_ForcedSphereObstacleOperator(FluidGridMPI*g, ArgumentParser&p, const Real*const u) : IF3D_SphereObstacleOperator(g,p,u)
  {
    _parseArguments(p);
  }

  void _parseArguments(ArgumentParser & parser) override;
  void update(const int stepID, const double t, const double dt, const Real *Uinf) override;
};

#endif /* defined(__IncompressibleFluids3D__IF3D_ForcedSphereObstacleOperator__) */
