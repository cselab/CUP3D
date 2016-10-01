//
//  CoordinatorPenalization.h
//  CubismUP_3D
//
//  Created by Christian Conti on 3/27/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef CubismUP_3D_CoordinatorPenalization_h
#define CubismUP_3D_CoordinatorPenalization_h

#include "GenericOperator.h"
#include "GenericCoordinator.h"
#include "IF3D_ObstacleVector.h"

class OperatorFadeOut : public GenericOperator
{
private:
	const Real extent[3];
	const int buffer;
public:
	OperatorFadeOut(const int buffer, const Real extent[3])
	: extent{extent[0],extent[1],extent[2]}, buffer(buffer) {}

	inline bool _is_touching(const BlockInfo& info, const Real h) const
	{
		Real max_pos[3],min_pos[3];
		info.pos(min_pos, 0, 0, 0);
		info.pos(max_pos, FluidBlock::sizeX-1, FluidBlock::sizeY-1, FluidBlock::sizeZ-1);
		// true: block within killing zone, false: block not within killing zone
		return  ( min_pos[0]<(0.+(2+buffer)*h ) ) ||
				( min_pos[1]<(0.+(2+buffer)*h ) ) ||
				( min_pos[2]<(0.+(2+buffer)*h ) ) ||
				( max_pos[0]>(extent[0]-(2+buffer)*h ) ) ||
				( max_pos[1]>(extent[1]-(2+buffer)*h ) ) ||
				( max_pos[2]>(extent[2]-(2+buffer)*h ) ) ;
	}

	void operator()(const BlockInfo& info, FluidBlock& block) const
	{
		const Real h = info.h_gridpoint;
		if(_is_touching(info,h))
		for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
		for(int iy=0; iy<FluidBlock::sizeY; ++iy)
		for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
			Real p[3];
			info.pos(p, ix, iy, iz);
			// >0 iff p+(1+buffer)*h > extent
			const Real argx1= max(0., ( p[0]-extent[0]+(2+buffer)*h) );
			const Real argy1= max(0., ( p[1]-extent[1]+(2+buffer)*h) );
			const Real argz1= max(0., ( p[2]-extent[2]+(2+buffer)*h) );
			// >0 iff (1+buffer)*h > p
			const Real argx2= max(0., ( 0.0 -p[0]     +(2+buffer)*h) );
			const Real argy2= max(0., ( 0.0 -p[1]     +(2+buffer)*h) );
			const Real argz2= max(0., ( 0.0 -p[2]     +(2+buffer)*h) );
			// max distance in killing zone 0 <= out <= (2+buffer)*h
			const Real out = max(max(max(argx1,argx2),max(argy1,argy2)),max(argz1,argz2));
			// 1 at buffer start, 0 at buffer end (2 grid points before border)
			const Real fade = max(Real(0.0), cos(0.5*M_PI* out/(buffer*h)));
			// smooth within killing zone (factor <= 1) and kill at very boundaries (factor < 0)
			b(ix,iy,iz).u = b(ix,iy,iz).u*fade;
			b(ix,iy,iz).v = b(ix,iy,iz).v*fade;
			b(ix,iy,iz).w = b(ix,iy,iz).w*fade;
			//b(ix,iy,iz).p = b(ix,iy,iz).p*fade;
		}
	}
};

class CoordinatorFadeOut : public GenericCoordinator
{
protected:
	const int buffer;
public:
    CoordinatorFadeOut(FluidGridMPI * grid, const int _buffer=5)
	: GenericCoordinator(grid), buffer(_buffer)
	{ }
	
	void operator()(const Real dt)
	{
		check((string)"FadeOut - start");

		const int N = vInfo.size();
		const Real h = grid->getH()
		const Real ext[3] = {
				h*grid->getBlocksPerDimension(0)*FluidBlock::sizeX,
				h*grid->getBlocksPerDimension(1)*FluidBlock::sizeY,
				h*grid->getBlocksPerDimension(2)*FluidBlock::sizeZ
		};
		#pragma omp parallel
		{
			OperatorFadeOut kernel(buffer,ext);
#pragma omp for schedule(static)
			for (int i=0; i<N; i++) {
				BlockInfo info = vInfo[i];
				FluidBlock& b = *(FluidBlock*)info.ptrBlock;
				kernel(info, b);
			}
		}

		check((string)"FadeOut - end");
	}
	
	string getName()
	{
		return "FadeOut";
	}
};

#endif
