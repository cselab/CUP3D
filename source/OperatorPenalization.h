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
	const Real uBody[2];
	const Real omegaBody;
	const Real centerOfMass[2];
	const double lambda;
	
public:
	OperatorPenalization(double dt, Real uSolid, Real vSolid, Real omegaBody, Real xCenterOfMass, Real yCenterOfMass, double lambda) : dt(dt), uBody{uSolid,vSolid}, omegaBody(omegaBody), centerOfMass{xCenterOfMass,yCenterOfMass}, lambda(lambda) {}
    ~OperatorPenalization() {}
    
    void operator()(const BlockInfo& info, FluidBlock& block) const
    {
		// this implementation considers that the Euler updates has already happened
		// do we need a finite state machine coordinating operators?
        for(int iy=0; iy<FluidBlock::sizeY; ++iy)
            for(int ix=0; ix<FluidBlock::sizeX; ++ix)
			{
				Real p[2] = {0,0};
                info.pos(p,ix,iy);
#ifndef _MOVING_FRAME_
                block(ix,iy).u = (block(ix,iy).u + dt * lambda * block(ix,iy).chi * (uBody[0] - omegaBody*(p[1]-centerOfMass[1]))) / (1 + dt * lambda * block(ix,iy).chi);
                block(ix,iy).v = (block(ix,iy).v + dt * lambda * block(ix,iy).chi * (uBody[1] + omegaBody*(p[0]-centerOfMass[0]))) / (1 + dt * lambda * block(ix,iy).chi);
#else
                block(ix,iy).u = (block(ix,iy).u + dt * lambda * block(ix,iy).chi * (uBody[0] - omegaBody*(p[1]-centerOfMass[1]))) / (1 + dt * lambda * block(ix,iy).chi);
                block(ix,iy).v = (block(ix,iy).v + dt * lambda * block(ix,iy).chi * (uBody[1] + omegaBody*(p[0]-centerOfMass[0]))) / (1 + dt * lambda * block(ix,iy).chi);
#endif
			}
    }
};


#endif
