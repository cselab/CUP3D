//
//  CoordinatorPressureGradient.h
//  CubismUP_3D
//
//  Created by Christian Conti on 5/16/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef __CubismUP_3D__CoordinatorPressureGradient__
#define __CubismUP_3D__CoordinatorPressureGradient__

#include "GenericCoordinator.h"
#include "OperatorPressureGradient.h"

class CoordinatorPressureGradient : public GenericCoordinator
{
protected:
	Real gradient[2];
	
public:
	CoordinatorPressureGradient(Real gradient[2], FluidGrid * grid) : GenericCoordinator(grid), gradient{gradient[0],gradient[1]}
	{
	}
	
	void operator()(const double dt)
	{
		check("gradient - start");
		
		BlockInfo * ary = &vInfo.front();
		const int N = vInfo.size();
		
#pragma omp parallel
		{
			OperatorPressureGradient kernel(gradient, dt);
			
#pragma omp for schedule(static)
			for (int i=0; i<N; i++)
				kernel(ary[i], *(FluidBlock*)ary[i].ptrBlock);
		}
		
		check("gradient - end");
	}
	
	string getName()
	{
		return "Pressure gradient";
	}
};

#endif /* defined(__CubismUP_3D__CoordinatorPressureGradient__) */
