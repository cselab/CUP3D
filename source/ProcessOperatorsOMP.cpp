//
//  ProcessOperatorsOMP.cpp
//  CubismUP_3D
//
//  Created by Christian Conti on 1/9/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#include "ProcessOperatorsOMP.h"
#include <cmath>

double findMaxUOMP(vector<BlockInfo>& myInfo, FluidGrid & grid)
{
	double maxU = 0;
	const int N = myInfo.size();
	
#pragma omp parallel for schedule(static) reduction(max:maxU)
	for(int i=0; i<N; i++)
	{
		BlockInfo info = myInfo[i];
		FluidBlock& b = *(FluidBlock*)info.ptrBlock;
		
		for(int iy=0; iy<FluidBlock::sizeY; ++iy)
			for(int ix=0; ix<FluidBlock::sizeX; ++ix)
			{
				maxU = max(maxU,(double)abs(b(ix,iy).u));
				maxU = max(maxU,(double)abs(b(ix,iy).v));
			}
	}
	
	return maxU;
};

void computeForcesFromVorticity(vector<BlockInfo>& myInfo, FluidGrid & grid, Real ub[2], Real oldAccVort[2], Real rhoS)
{
	Real mU = 0;
	Real mV = 0;
	Real mass = 0;
	const int N = myInfo.size();
	
#pragma omp parallel for schedule(static) reduction(+:mU) reduction(+:mV) reduction(+:mass)
	for(int i=0; i<N; i++)
	{
		BlockInfo info = myInfo[i];
		FluidBlock& b = *(FluidBlock*)info.ptrBlock;
		
		Real h = info.h_gridpoint;
		
		for(int iy=0; iy<FluidBlock::sizeY; ++iy)
			for(int ix=0; ix<FluidBlock::sizeX; ++ix)
			{
				Real p[2];
				info.pos(p, ix, iy);
				
				mU += - (1-b(ix,iy).chi) * p[1] * b(ix,iy).tmp * b(ix,iy).rho;
				mV -= - (1-b(ix,iy).chi) * p[0] * b(ix,iy).tmp * b(ix,iy).rho;
				mass += b(ix,iy).chi * rhoS;
			}
	}
	
	ub[0] += (mU-oldAccVort[0]) / mass;
	ub[1] += (mV-oldAccVort[1]) / mass;
	oldAccVort[0] = mU;
	oldAccVort[1] = mV;
}
