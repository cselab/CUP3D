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
	const Real centerOfMass[3];
	const double lambda;
	
public:
	OperatorPenalization(double dt, Real uSolid, Real vSolid, Real wSolid, Real xCenterOfMass, Real yCenterOfMass, Real zCenterOfMass, double lambda) : dt(dt), uBody{uSolid,vSolid,wSolid}, centerOfMass{xCenterOfMass,yCenterOfMass,zCenterOfMass}, lambda(lambda) {}
    ~OperatorPenalization() {}
    
    void operator()(const BlockInfo& info, FluidBlock& block) const
	{		
		// this implementation considers that the Euler updates has already happened
		// do we need a finite state machine coordinating operators?
		for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
		for(int iy=0; iy<FluidBlock::sizeY; ++iy)
            for(int ix=0; ix<FluidBlock::sizeX; ++ix)
			{
				block(ix,iy,iz).u = (block(ix,iy,iz).u + dt * lambda * block(ix,iy,iz).chi * uBody[0]) / (1 + dt * lambda * block(ix,iy,iz).chi);
				block(ix,iy,iz).v = (block(ix,iy,iz).v + dt * lambda * block(ix,iy,iz).chi * uBody[1]) / (1 + dt * lambda * block(ix,iy,iz).chi);
				block(ix,iy,iz).w = (block(ix,iy,iz).w + dt * lambda * block(ix,iy,iz).chi * uBody[2]) / (1 + dt * lambda * block(ix,iy,iz).chi);
			}
    }
};


#endif
