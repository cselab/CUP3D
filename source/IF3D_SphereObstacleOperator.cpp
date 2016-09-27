/*
 *  IF3D_SphereObstacleOperator.cpp
 *  IncompressibleFluids3D
 *
 *  Created by Diego Rossinelli on 8/6/10.
 *  Copyright 2010 ETH Zurich. All rights reserved.
 *
 */

#include <limits>
#include <cstring>

#include "IF3D_SphereObstacleOperator.h"
#include "IF3D_ObstacleLibrary.h"

void IF3D_SphereObstacleOperator::create(const int step_id,const Real time, const Real dt, const Real *Uinf)
{
    for(auto & entry : obstacleBlocks)
        delete entry.second;
    obstacleBlocks.clear();

    SphereObstacle::FillBlocks kernel(radius,vInfo[0].h_gridpoint,position);
    for(int i=0; i<vInfo.size(); i++) {
    	BlockInfo info = vInfo[i];
        const auto pos = obstacleBlocks.find(info.blockID);
		if(kernel._is_touching(info)) { //position of sphere + radius + 2*h safety
			assert(obstacleBlocks.find(info.blockID) == obstacleBlocks.end());
			obstacleBlocks[info.blockID] = new ObstacleBlock;
			obstacleBlocks[info.blockID]->clear(); //memset 0
		}
    }

	const int nthreads = omp_get_max_threads();
	vector<surfaceBlocks> dataPerThread(nthreads);

#pragma omp parallel
	{
		SphereObstacle::FillBlocks fill(radius,vInfo[0].h_gridpoint,position);
		const int tid = omp_get_thread_num();

#pragma omp for schedule(static)
		for(int i=0; i<vInfo.size(); i++) {
			BlockInfo info = vInfo[i];
			auto pos = obstacleBlocks.find(info.blockID);
			if(pos == obstacleBlocks.end()) continue;
			FluidBlock& b = *(FluidBlock*)info.ptrBlock;
			fill(info, b, pos->second, &dataPerThread[tid]);
		}
	}
	surfData.finalizeOnGrid(dataPerThread);
}

void IF3D_SphereObstacleOperator::_parseArguments(ArgumentParser & parser)
{
	//obstacleop parses x,y,z,quats and length!
	IF3D_ObstacleOperator::_parseArguments(parser);
    radius = .5*length;
}
