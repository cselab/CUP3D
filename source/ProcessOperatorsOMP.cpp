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
		
		for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
		for(int iy=0; iy<FluidBlock::sizeY; ++iy)
			for(int ix=0; ix<FluidBlock::sizeX; ++ix)
			{
				maxU = max(maxU,(double)abs(b(ix,iy,iz).u));
				maxU = max(maxU,(double)abs(b(ix,iy,iz).v));
				maxU = max(maxU,(double)abs(b(ix,iy,iz).w));
			}
	}
	
	return maxU;
};

void computeForcesFromVorticity(vector<BlockInfo>& myInfo, FluidGrid & grid, Real ub[3], Real oldAccVort[3], Real rhoS)
{
	Real mU = 0;
	Real mV = 0;
	Real mW = 0;
	Real mass = 0;
	const int N = myInfo.size();
	
#pragma omp parallel for schedule(static) reduction(+:mU) reduction(+:mV) reduction(+:mW) reduction(+:mass)
	for(int i=0; i<N; i++)
	{
		BlockInfo info = myInfo[i];
		FluidBlock& b = *(FluidBlock*)info.ptrBlock;
		
		Real h = info.h_gridpoint;
		
		for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
		for(int iy=0; iy<FluidBlock::sizeY; ++iy)
			for(int ix=0; ix<FluidBlock::sizeX; ++ix)
			{
				Real p[3];
				info.pos(p, ix, iy, iz);
				
				mU += - (1-b(ix,iy,iz).chi) * p[1] * b(ix,iy,iz).tmp * b(ix,iy,iz).rho;
				mV -= - (1-b(ix,iy,iz).chi) * p[0] * b(ix,iy,iz).tmp * b(ix,iy,iz).rho;
				mW  = 0;
				mass += b(ix,iy,iz).chi * rhoS;
				abort(); // something is missing here for the vorticity
			}
	}
	
	ub[0] += (mU-oldAccVort[0]) / mass;
	ub[1] += (mV-oldAccVort[1]) / mass;
	ub[2] += (mV-oldAccVort[2]) / mass;
	oldAccVort[0] = mU;
	oldAccVort[1] = mV;
	oldAccVort[2] = mW;
}
