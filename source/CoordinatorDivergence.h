//
//  CoordinatorDivergence.h
//  CubismUP_3D
//
//  Created by Christian Conti on 3/15/16.
//  Copyright © 2016 ETHZ. All rights reserved.
//

#ifndef CoordinatorDivergence_h
#define CoordinatorDivergence_h

#include "GenericCoordinator.h"
#include "OperatorDivergence.h"

template <typename Lab>
class CoordinatorDivergence : public GenericCoordinator
{
protected:
	
public:
	CoordinatorDivergence(FluidGridMPI * grid) : GenericCoordinator(grid)
	{
	}
	
	void operator()(const double dt)
	{
		check("divergence - start");
		
		OperatorDivergence kernel(1);
		compute(kernel);
		
		check("divergence - end");
	}
	
	string getName()
	{
		return "Divergence";
	}
};


#endif /* CoordinatorDivergence_h */
