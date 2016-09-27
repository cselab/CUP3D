/*
 *  IF3D_SphereObstacleOperator.h
 *  IncompressibleFluids3D
 *
 *  Created by Diego Rossinelli on 8/6/10.
 *  Copyright 2010 ETH Zurich. All rights reserved.
 *
 *	compute the characteristic function of an obstacle and
 *	store it into tmp. Tested!
 *
 *	characteristic_function():
 *	IN: 
 *	OUT: tmp
 *
 */
#pragma once

#include "IF3D_ObstacleOperator.h"

#include <cmath>

class IF3D_SphereObstacleOperator: public IF3D_ObstacleOperator
{
	Real radius;
	
public:
	
 IF3D_SphereObstacleOperator(FluidGridMPI * grid, ArgumentParser & parser) //const Real radius, const double position[3], const Real smoothing_length=-1):
	: IF3D_ObstacleOperator(grid, parser)//, radius(radius), smoothing_length(smoothing_length)
	{
	 	 _parseArguments(parser);
	}

 	void create(const int step_id,const Real time, const Real dt, const Real *Uinf) override;
    void _parseArguments(ArgumentParser & parser) override;
};
