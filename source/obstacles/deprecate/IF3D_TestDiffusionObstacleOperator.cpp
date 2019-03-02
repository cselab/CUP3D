//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by SV (vermas@fau.edu).
//

#include "obstacles/IF3D_TestDiffusionObstacleOperator.h"
#include "obstacles/IF3D_ObstacleLibrary.h"

#include "Cubism/ArgumentParser.h"

IF3D_TestDiffusionObstacleOperator::IF3D_TestDiffusionObstacleOperator(
    SimulationData& s, ArgumentParser &p )
    : IF3D_ObstacleOperator(s, p), sigma(p("-sigma").asDouble(0.0))
{
	printf("Created IF3D_TestDiffusionObstacleOperator with sigma %f\n", sigma);
}


void IF3D_TestDiffusionObstacleOperator::create(const int step_id,const double time, const double dt, const Real *Uinf)
{

	if(step_id>0 || time>0.0) return;

	//this writes the initial u field
  #pragma omp parallel
	{
    #pragma omp for schedule(dynamic)
		for(size_t i=0; i<vInfo.size(); i++) {

			BlockInfo info = vInfo[i];
			FluidBlock& block = *(FluidBlock*)info.ptrBlock;

			for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
				for(int iy=0; iy<FluidBlock::sizeY; ++iy)
					for(int ix=0; ix<FluidBlock::sizeX; ++ix) {

						block(ix,iy,iz).u = 0;//.05*cos(x+0.00*M_PI)*cos(y+0.75*M_PI)*cos(z+1.50*M_PI);
						block(ix,iy,iz).v = 0;//.05*cos(x+0.25*M_PI)*cos(y+1.00*M_PI)*cos(z+1.75*M_PI);
						block(ix,iy,iz).w = 0;//.05*cos(x+0.50*M_PI)*cos(y+1.25*M_PI)*cos(z+2.00*M_PI);
					}
		}
	}
}
