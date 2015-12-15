//
//  OperatorPenalization.h
//  CubismUP_3D
//
//	Operates on
//		u, v
//
//  Created by Christian Conti on 1/7/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef CubismUP_3D_OperatorPenalization_h
#define CubismUP_3D_OperatorPenalization_h

#include "GenericOperator.h"

class OperatorPenalization : public GenericOperator
{
private:
	const double dt;
	const Real uBody[3];
	Shape * shape;
	const double lambda;
	
public:
	OperatorPenalization(double dt, Real uSolid, Real vSolid, Real wSolid, Shape * shape, double lambda) : dt(dt), uBody{uSolid,vSolid,wSolid}, shape(shape), lambda(lambda) {}
	~OperatorPenalization() {}
	
	void operator()(const BlockInfo& info, FluidBlock& block) const
	{
		// this implementation considers that the Euler updates has already happened
		// do we need a finite state machine coordinating operators?
		for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
			for(int iy=0; iy<FluidBlock::sizeY; ++iy)
				for(int ix=0; ix<FluidBlock::sizeX; ++ix)
				{
					Real dtChiLambda = dt * lambda * block(ix,iy,iz).chi;
					block(ix,iy,iz).u = (block(ix,iy,iz).u + dtChiLambda * uBody[0]) / (1. + dtChiLambda);
					block(ix,iy,iz).v = (block(ix,iy,iz).v + dtChiLambda * uBody[1]) / (1. + dtChiLambda);
					block(ix,iy,iz).w = (block(ix,iy,iz).w + dtChiLambda * uBody[2]) / (1. + dtChiLambda);
					
					Real p[3];
					info.pos(p, ix, iy, iz);
					
					Real chi = shape->chi(p, info.h_gridpoint);
					block(ix,iy,iz).chi = chi;
					block(ix,iy,iz).rho = shape->rho(p, info.h_gridpoint, chi);
				}
	}
};


#endif
