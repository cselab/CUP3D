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
