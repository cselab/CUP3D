//
//  IF3D_ForcedSphereObstacleOperator.cpp
//  IncompressibleFluids3D
//
//  Created by Wim van Rees on 8/20/13.
//
//

#include <limits>
#include "IF3D_ForcedSphereObstacleOperator.h"
#include "IF3D_ObstacleLibrary.h"


void IF3D_ForcedSphereObstacleOperator::_parseArguments(ArgumentParser & parser)
{
	IF3D_SphereObstacleOperator::_parseArguments(parser);

  parser.set_strict_mode();
  const Real xvel = parser("-xvel").asDouble();
  parser.unset_strict_mode();
  const Real yvel = parser("-yvel").asDouble(0.);
  const Real zvel = parser("-zvel").asDouble(0.);
  accel_decel = parser("-accel").asBool(false);
  umax = xvel;
  tmax = parser("-T").asDouble(1.);

	this->transVel[0] = xvel;
	this->transVel[1] = yvel;
	this->transVel[2] = zvel;
}

void IF3D_ForcedSphereObstacleOperator::update(const int stepID, const Real t, const Real dt, const Real *Uinf)
{
  if(accel_decel)
  {
    if(t<tmax)
      this->transVel[0] = umax*t/tmax;
    else if (t<2*tmax)
      this->transVel[0] = umax*(2*tmax-t)/tmax;
    else
      this->transVel[0] = 0;
  }
	// update position and angles
	IF3D_ObstacleOperator::update(stepID,t, dt, Uinf);
}
