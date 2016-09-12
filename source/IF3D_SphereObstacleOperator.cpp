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

void IF3D_SphereObstacleOperator::create()
{
	const int nthreads = omp_get_max_threads();
	vector<surfaceBlocks> dataPerThread(nthreads);

#pragma omp parallel
	{
		SphereObstacle::FillBlocks fill(radius,vInfo[0].h_gridpoint,position);
		const int tid = omp_get_thread_num();

#pragma omp for schedule(static)
		for(int i=0; i<vInfo.size(); i++) {
			auto pos = obstacleBlocks.find(i);
			if(pos == obstacleBlocks.end()) continue;
			FluidBlock& b = *(FluidBlock*)ary[i].ptrBlock;

			lab.load(ary[i], 0);
			double tmp[4];
			finalize(lab, ary[i], b, pos->second, dataPerThread(tid));
		}
	}
	surfData.finalizeOnGrid(dataPerThread);
}

void IF3D_SphereObstacleOperator::_parseArguments(ArgumentParser & parser)
{
	IF2D_ObstacleOperator::_parseArguments(parser);

    parser.set_strict_mode();
    length = parser("-D").asDouble();
    radius = .5*length;
}
