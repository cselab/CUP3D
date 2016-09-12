//
//  IF3D_ForcedSphereObstacleOperator.cpp
//  IncompressibleFluids3D
//
//  Created by Wim van Rees on 8/20/13.
//
//

#include <limits>
#include "IF3D_ForcedSphereObstacleOperator.h"


void IF3D_ForcedSphereObstacleOperator::_parseArguments(ArgumentParser & parser)
{
	IF3D_SphereObstacleOperator::_parseArguments(parser);

    parser.set_strict_mode();
    Real xvel = parser("-xvel").asDouble();
    parser.unset_strict_mode();
    Real yvel = parser("-yvel").asDouble(0.);
    Real zvel = parser("-zvel").asDouble(0.);

	this->transVel[0] = xvel;
	this->transVel[1] = yvel;
	this->transVel[2] = zvel;
}
