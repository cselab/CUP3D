//
//  OperatorComputeShape.h
//  CubismUP_3D
//
//	Operates on
//		chi, rho
//
//  Created by Christian Conti on 1/22/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef CubismUP_3D_OperatorComputeShape_h
#define CubismUP_3D_OperatorComputeShape_h

#include "GenericOperator.h"

class OperatorComputeShape : public GenericOperator
{
private:
	Shape * shape;
	
public:
	OperatorComputeShape(Shape * shape) : shape(shape) {}
	~OperatorComputeShape() {}
	
	void operator()(const BlockInfo& info, FluidBlock& block) const
	{
		for(int iy=0; iy<FluidBlock::sizeY; ++iy)
			for(int ix=0; ix<FluidBlock::sizeX; ++ix)
			{
				Real p[2];
				info.pos(p, ix, iy);
                
				Real chi = shape->chi(p, info.h_gridpoint);
				block(ix,iy).chi = chi;
				block(ix,iy).rho = shape->getRhoS()*chi + 1.*(1.-chi);//shape->rho(p, info.h_gridpoint);
			}
	}
};

#endif
