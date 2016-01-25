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
	CoordinatorPenalizationFixed(Shape * shape, Real * lambda, FluidGridMPI * grid) : GenericCoordinator(grid), shape(shape), lambda(lambda)
	{
	}
	
	void operator()(const double dt)
	{
		check("penalization - start");
		
		BlockInfo * ary = &vInfo.front();
		const int N = vInfo.size();
		
#pragma omp parallel
		{
			OperatorPenalization kernel(dt, 0, 0, 0, shape, *lambda);
			
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
	Real *uBody, *vBody, *wBody;
	Shape * shape;
	Real * lambda;
	
public:
	CoordinatorPenalization(Real * uBody, Real * vBody, Real * wBody, Shape * shape, Real * lambda, FluidGridMPI * grid) : GenericCoordinator(grid), uBody(uBody), vBody(vBody), wBody(wBody), shape(shape), lambda(lambda)
	{
	}
	
	void operator()(const double dt)
	{
		check("penalization - start");
		
		BlockInfo * ary = &vInfo.front();
		const int N = vInfo.size();
		
#pragma omp parallel
		{
			OperatorPenalization kernel(dt, *uBody, *vBody, *wBody, shape, *lambda);
			
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
