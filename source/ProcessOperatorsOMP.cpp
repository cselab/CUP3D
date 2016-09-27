//
//  ProcessOperatorsOMP.cpp
//  CubismUP_3D
//
//  Created by Christian Conti on 1/9/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#include "ProcessOperatorsOMP.h"
#include <cmath>

Real findMaxUOMP(vector<BlockInfo>& myInfo, FluidGridMPI & grid, const Real* const uInf)
{
	Real maxU = 0;
	const int N = myInfo.size();
	
#pragma omp parallel for schedule(static) reduction(max:maxU)
	for(int i=0; i<N; i++)
	{
		BlockInfo info = myInfo[i];
		FluidBlock& b = *(FluidBlock*)info.ptrBlock;
		
		for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
		for(int iy=0; iy<FluidBlock::sizeY; ++iy)
		for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
			maxU = max(maxU,(Real)abs(b(ix,iy,iz).u +uInf[0]));
			maxU = max(maxU,(Real)abs(b(ix,iy,iz).v +uInf[1]));
			maxU = max(maxU,(Real)abs(b(ix,iy,iz).w +uInf[2]));
		}
	}
	
	return maxU;
};
