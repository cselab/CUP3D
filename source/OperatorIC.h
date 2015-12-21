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
		const Real dh = info.h_gridpoint;
		for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
			for(int iy=0; iy<FluidBlock::sizeY; ++iy)
				for(int ix=0; ix<FluidBlock::sizeX; ++ix)
				{
					Real p[3];
					info.pos(p, ix, iy, iz);
					
					block(ix,iy,iz).u = uinf;
					block(ix,iy,iz).v = 0;
					block(ix,iy,iz).w = 0;
					block(ix,iy,iz).chi = shape->chi(p, dh);
					
					// assume fluid with density 1
					block(ix,iy,iz).rho = shape->rho(p, dh, block(ix,iy,iz).chi);
					
					block(ix,iy,iz).p = 0;
					block(ix,iy,iz).divU = 0;
					block(ix,iy,iz).pOld = 0;
					
					block(ix,iy,iz).tmpU = 0;
					block(ix,iy,iz).tmpV = 0;
					block(ix,iy,iz).tmpW = 0;
					block(ix,iy,iz).tmp  = 0;
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
		
		for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
			for(int iy=0; iy<FluidBlock::sizeY; ++iy)
				for(int ix=0; ix<FluidBlock::sizeX; ++ix)
				{
					Real p[3];
					info.pos(p, ix, iy, iz);
					
					Real x = p[0] - d*.5;
					Real y = p[1] - d*2.5;
					Real z = p[2] - d*.5;
					/*
					 Real eta = -.1*a*cos(2*M_PI*x/d);
					 block(ix,iy).rho = 2. + tanh((y-eta)/(.01*d));
					 */
					Real eta = -a*.25*cos(2*M_PI*x/d)*cos(2*M_PI*z/d);
					//block(ix,iy).rho = (1.+rhoS)/2. + ((1.+rhoS)/2.-1.)*tanh((y-eta)/(.01*d));
					block(ix,iy,iz).rho = (1.+rhoS)/2. + ((1.+rhoS)/2.-1.)*tanh((y-eta)/info.h_gridpoint);
					block(ix,iy,iz).u = 0;
					block(ix,iy,iz).v = 0;
					block(ix,iy,iz).w = 0;
					block(ix,iy,iz).chi = .5+.5*tanh((y-eta)/info.h_gridpoint);
					
					block(ix,iy,iz).p = 0;
					block(ix,iy,iz).divU = 0;
					block(ix,iy,iz).pOld = 0;
					
					block(ix,iy,iz).tmpU = 0;
					block(ix,iy,iz).tmpV = 0;
					block(ix,iy,iz).tmpW = 0;
					block(ix,iy,iz).tmp  = 0;
				}
	}
};


#endif
