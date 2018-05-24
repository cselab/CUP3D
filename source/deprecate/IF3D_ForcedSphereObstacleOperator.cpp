//
//  CubismUP_3D
//
//  Written by Guido Novati ( novatig@ethz.ch ).
//  This file started as an extension of code written by Wim van Rees
//  Copyright (c) 2017 ETHZ. All rights reserved.
//

#include <limits>
#include "IF3D_ForcedSphereObstacleOperator.h"
#include "IF3D_ObstacleLibrary.h"


void IF3D_ForcedSphereObstacleOperator::_parseArguments(ArgumentParser & parser)
{
  IF3D_SphereObstacleOperator::_parseArguments(parser);
  bForcedInSimFrame[0] = true;
  bForcedInSimFrame[1] = true;
  bForcedInSimFrame[2] = true;
}
