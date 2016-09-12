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

class OperatorPressureGradient : public GenericOperator
{
private:
	double dt;
	Real gradient[3];

public:
	OperatorPressureGradient(Real gradient[3], double dt) : dt(dt), gradient{gradient[0],gradient[1],gradient[2]}
	{
	}

	~OperatorPressureGradient() {}

	void operator()(const BlockInfo& info, FluidBlock& block) const
	{
		for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
			for(int iy=0; iy<FluidBlock::sizeY; ++iy)
			for(int ix=0; ix<FluidBlock::sizeX; ++ix)
			{
				// assumes [0,1[^2 domain

				double p[3];
				info.pos(p,ix,iy,iz);

				block(ix,iy,iz).u += dt*gradient[0]*(1-p[0]);
				block(ix,iy,iz).v += dt*gradient[1]*(1-p[1]);
				block(ix,iy,iz).w += dt*gradient[2]*(1-p[2]);
			}
	}
};

class CoordinatorPressureGradient : public GenericCoordinator
{
protected:
	Real gradient[3];
	
public:
	CoordinatorPressureGradient(Real gradient[3], FluidGridMPI * grid) : GenericCoordinator(grid), gradient{gradient[0],gradient[1],gradient[2]}
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
