//
//  CoordinatorComputeShape.h
//  CubismUP_3D
//
//  Created by Christian Conti on 3/30/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef CubismUP_3D_CoordinatorComputeShape_h
#define CubismUP_3D_CoordinatorComputeShape_h

#include "GenericCoordinator.h"
#include "OperatorComputeShape.h"
#include "Shape.h"

class CoordinatorComputeShape : public GenericCoordinator
{
protected:
	Real *uBody, *vBody, *omegaBody;
	Shape * shape;
    
public:
	CoordinatorComputeShape(Real * uBody, Real * vBody, Real * omegaBody, Shape * shape, FluidGrid * grid) : GenericCoordinator(grid), uBody(uBody), vBody(vBody), omegaBody(omegaBody), shape(shape)
	{
	}
	
	void operator()(const double dt)
	{
		check("shape - start");
		
		BlockInfo * ary = &vInfo.front();
		const int N = vInfo.size();
		
		Real ub[2] = { *uBody, *vBody };
		shape->updatePosition(ub, *omegaBody, dt);
		
		Real domainSize[2] = { grid->getBlocksPerDimension(0)*FluidBlock::sizeX*vInfo[0].h_gridpoint,
							   grid->getBlocksPerDimension(1)*FluidBlock::sizeY*vInfo[0].h_gridpoint};
		Real p[2] = {0,0};
		shape->getPosition(p);
		
		if (p[0]<0 || p[0]>domainSize[0] || p[1]<0 || p[1]>domainSize[1])
			exit(0);
		
#pragma omp parallel
		{
			OperatorComputeShape kernel(shape);
			
#pragma omp for schedule(static)
			for(int i=0; i<N; i++)
				kernel(ary[i], *(FluidBlock*)ary[i].ptrBlock);
		}
		
		check("shape - end");
	}
	
	string getName()
	{
		return "ComputeShape";
	}
};

#endif
