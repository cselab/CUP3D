//
//  OperatorIC.h
//  CubismUP_3D
//
//	Operates on
//		chi, u, v, rho, p, pOld
//
//  Created by Christian Conti on 1/7/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef CubismUP_3D_OperatorIC_h
#define CubismUP_3D_OperatorIC_h

#include "GenericOperator.h"

class OperatorIC : public GenericOperator
{
private:
	Shape * shape;
	const double uinf;
	
public:
	OperatorIC(Shape * shape, const double uinf) : shape(shape), uinf(uinf) {}
	
	~OperatorIC() {}
	
	void operator()(const BlockInfo& info, FluidBlock& block) const
	{
		for(int iy=0; iy<FluidBlock::sizeY; ++iy)
			for(int ix=0; ix<FluidBlock::sizeX; ++ix)
			{
				Real p[2];
				info.pos(p, ix, iy);
				
				block(ix,iy).u = uinf;
				block(ix,iy).v = 0;//uinf/1000.*sin(p[0]*100);
				block(ix,iy).chi = shape->chi(p, info.h_gridpoint);
				
				// assume fluid with density 1
				block(ix,iy).rho = shape->rho(p, info.h_gridpoint);
				
				block(ix,iy).p = 0;
				block(ix,iy).divU = 0;
				block(ix,iy).pOld = 0;
				
				block(ix,iy).tmpU = 0;
				block(ix,iy).tmpV = 0;
				block(ix,iy).tmp  = 0;
			}
	}
};

class OperatorIC_RT : public GenericOperator
{
protected:
	const double rhoS;
	
public:
	OperatorIC_RT(const double rhoS) : rhoS(rhoS) {}
	
	~OperatorIC_RT() {}
	
	void operator()(const BlockInfo& info, FluidBlock& block) const
	{
		Real d = .25;
		Real a = .05;//.25
		
		for(int iy=0; iy<FluidBlock::sizeY; ++iy)
			for(int ix=0; ix<FluidBlock::sizeX; ++ix)
			{
				Real p[2];
				info.pos(p, ix, iy);
				
				Real y = p[1] - d*2.;
				Real x = p[0] - d*.5;
				/*
				Real eta = -.1*a*cos(2*M_PI*x/d);
				block(ix,iy).rho = 2. + tanh((y-eta)/(.01*d));
				 */
				Real eta = -a*.25*cos(2*M_PI*x/d);
				//block(ix,iy).rho = (1.+rhoS)/2. + ((1.+rhoS)/2.-1.)*tanh((y-eta)/(.01*d));
				block(ix,iy).rho = (1.+rhoS)/2. + ((1.+rhoS)/2.-1.)*tanh((y-eta)/info.h_gridpoint);
				block(ix,iy).u = 0;
				block(ix,iy).v = 0;
				block(ix,iy).chi = .5+.5*tanh((y-eta)/info.h_gridpoint);
				
				block(ix,iy).p = 0;
				block(ix,iy).divU = 0;
				block(ix,iy).pOld = 0;
				
				block(ix,iy).tmpU = 0;
				block(ix,iy).tmpV = 0;
				block(ix,iy).tmp  = 0;
			}
	}
};


#endif
