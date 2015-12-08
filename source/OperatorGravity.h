//
//  OperatorGravity.h
//  CubismUP_3D
//
//	Operates on
//		u, v
//
//  Created by Christian Conti on 1/22/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef CubismUP_3D_OperatorGravity_h
#define CubismUP_3D_OperatorGravity_h

#include "GenericOperator.h"

class OperatorGravity : public GenericOperator
{
private:
	const double dt;
	const Real g[3];
	
public:
	OperatorGravity(Real g[3], double dt) : dt(dt), g{g[0],g[1],g[2]} {}
	~OperatorGravity() {}
	
	void operator()(const BlockInfo& info, FluidBlock& block) const
	{
		for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
			for(int iy=0; iy<FluidBlock::sizeY; ++iy)
			for(int ix=0; ix<FluidBlock::sizeX; ++ix)
			{
				block(ix,iy,iz).u += dt * g[0];
				block(ix,iy,iz).v += dt * g[1];
				block(ix,iy,iz).w += dt * g[2];
			}
	}
};

#endif
