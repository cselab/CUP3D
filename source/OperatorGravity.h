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
	const Real g[2];
	
public:
	OperatorGravity(Real g[2], double dt) : dt(dt), g{g[0],g[1]} {}
	~OperatorGravity() {}
	
	void operator()(const BlockInfo& info, FluidBlock& block) const
	{
		for(int iy=0; iy<FluidBlock::sizeY; ++iy)
			for(int ix=0; ix<FluidBlock::sizeX; ++ix)
			{
				Real p[2];
				info.pos(p, ix, iy);
				
				block(ix,iy).u += dt * g[0];
				block(ix,iy).v += dt * g[1];
			}
	}
};

#endif
