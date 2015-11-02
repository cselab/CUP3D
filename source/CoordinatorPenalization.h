//
//  CoordinatorPenalization.h
//  CubismUP_3D
//
//  Created by Christian Conti on 3/27/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef CubismUP_3D_CoordinatorPenalization_h
#define CubismUP_3D_CoordinatorPenalization_h

#include "GenericCoordinator.h"
#include "OperatorPenalization.h"
#include "Shape.h"

class CoordinatorPenalizationFixed : public GenericCoordinator
{
protected:
	Shape * shape;
	Real * lambda;
	
public:
	CoordinatorPenalizationFixed(Shape * shape, Real * lambda, FluidGrid * grid) : GenericCoordinator(grid), shape(shape), lambda(lambda)
	{
	}
	
	void operator()(const double dt)
	{
		check("penalization - start");
		
		BlockInfo * ary = &vInfo.front();
		const int N = vInfo.size();
		
#pragma omp parallel
		{
			Real com[3];
			shape->getPosition(com);
			OperatorPenalization kernel(dt, 0, 0, 0, com[0], com[1], com[2], *lambda);
			
#pragma omp for schedule(static)
			for(int i=0; i<N; i++)
				kernel(ary[i], *(FluidBlock*)ary[i].ptrBlock);
		}
		
		check("penalization - end");
	}
	
	string getName()
	{
		return "Penalization";
	}
};

class CoordinatorPenalization : public GenericCoordinator
{
protected:
	Real *uBody, *vBody, *wBody, *omegaBody;
	Shape * shape;
	Real * lambda;
	
public:
	CoordinatorPenalization(Real * uBody, Real * vBody, Real * wBody, Real * omegaBody, Shape * shape, Real * lambda, FluidGrid * grid) : GenericCoordinator(grid), uBody(uBody), vBody(vBody), wBody(wBody), omegaBody(omegaBody), shape(shape), lambda(lambda)
	{
	}
	
	void operator()(const double dt)
	{
		check("penalization - start");
		
		BlockInfo * ary = &vInfo.front();
		const int N = vInfo.size();
		
#pragma omp parallel
		{
			Real com[3];
			shape->getPosition(com);
			OperatorPenalization kernel(dt, *uBody, *vBody, *wBody, *omegaBody, com[0], com[1], *lambda);
			
#pragma omp for schedule(static)
			for(int i=0; i<N; i++)
				kernel(ary[i], *(FluidBlock*)ary[i].ptrBlock);
		}
		
		check("penalization - end");
	}
	
	string getName()
	{
		return "Penalization";
	}
};

#endif
