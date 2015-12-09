//
//  CoordinatorVorticity.h
//  CubismUP_3D
//
//  Created by Christian Conti on 12/8/15.
//  Copyright © 2015 ETHZ. All rights reserved.
//

#ifndef CubismUP_3D_CoordinatorVorticity_h
#define CubismUP_3D_CoordinatorVorticity_h

#include "GenericCoordinator.h"
#include "OperatorVorticity.h"

template <typename Lab>
class CoordinatorVorticity : public GenericCoordinator
{
protected:
	
	inline void reset()
	{
		const int N = vInfo.size();
		
#pragma omp parallel for schedule(static)
		for(int i=0; i<N; i++)
		{
			BlockInfo info = vInfo[i];
			FluidBlock& b = *(FluidBlock*)info.ptrBlock;
			
			for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
			for(int iy=0; iy<FluidBlock::sizeY; ++iy)
			for(int ix=0; ix<FluidBlock::sizeX; ++ix)
			{
				b(ix,iy,iz).tmpU = 0;
				b(ix,iy,iz).tmpV = 0;
				b(ix,iy,iz).tmpW = 0;
				b(ix,iy,iz).tmp = 0;
			}
		}
	};
	
	inline void vorticity()
	{
		BlockInfo * ary = &vInfo.front();
		const int N = vInfo.size();
		
#pragma omp parallel
		{
			OperatorVorticity kernel;
			
			Lab mylab;
			mylab.prepare(*grid, kernel.stencil_start, kernel.stencil_end, false);
			
#pragma omp for schedule(static)
			for (int i=0; i<N; i++)
			{
				mylab.load(ary[i], 0);
				kernel(mylab, ary[i], *(FluidBlock*)ary[i].ptrBlock);
			}
		}
	}
	
public:
	CoordinatorVorticity(FluidGrid * grid) : GenericCoordinator(grid)
	{
	}
	
	void operator()(const double dt)
	{
		check("vorticity - start");
		
		reset();
		vorticity();
		
		check("vorticity - end");
	}
	
	string getName()
	{
		return "Vorticity";
	}
};


#endif /* CoordinatorVorticity_h */
