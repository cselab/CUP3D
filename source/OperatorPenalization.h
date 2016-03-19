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
		Real com[3];
		shape->getCenterOfMass(com);
		
		Real omegaBody[3];
		shape->getAngularVelocity(omegaBody);
		
		// this implementation considers that the Euler updates has already happened
		// do we need a finite state machine coordinating operators?
		for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
			for(int iy=0; iy<FluidBlock::sizeY; ++iy)
				for(int ix=0; ix<FluidBlock::sizeX; ++ix)
				{
					Real p[3], r[3];
					info.pos(p,ix,iy,iz);
					r[0] = p[0] - com[0];
					r[1] = p[1] - com[1];
					r[2] = p[2] - com[2];
					
					Real urBody[3] = {	omegaBody[1]*r[2]-omegaBody[2]*r[1],
										omegaBody[2]*r[0]-omegaBody[0]*r[2],
										omegaBody[0]*r[1]-omegaBody[1]*r[0] };
					
					Real dtChiLambda = dt * lambda * block(ix,iy,iz).chi;
					
//#ifdef _DLM_
//					block(ix,iy,iz).u += dtChiLambda * (uBody[0]+urBody[0] - block(ix,iy,iz).u);
//					block(ix,iy,iz).v += dtChiLambda * (uBody[1]+urBody[1] - block(ix,iy,iz).v);
//					block(ix,iy,iz).w += dtChiLambda * (uBody[2]+urBody[2] - block(ix,iy,iz).w);
//#else
					block(ix,iy,iz).u = (block(ix,iy,iz).u + dtChiLambda * (uBody[0]+urBody[0])) / (1. + dtChiLambda);
					block(ix,iy,iz).v = (block(ix,iy,iz).v + dtChiLambda * (uBody[1]+urBody[1])) / (1. + dtChiLambda);
					block(ix,iy,iz).w = (block(ix,iy,iz).w + dtChiLambda * (uBody[2]+urBody[2])) / (1. + dtChiLambda);
//#endif
					
					Real chi = shape->chi(p, info.h_gridpoint);
					block(ix,iy,iz).chi = chi;
					block(ix,iy,iz).rho = shape->rho(p, info.h_gridpoint, chi);
				}
	}
};


class OperatorPenalizationMovingFrame : public GenericOperator
{
private:
	const double dt;
	const Real uBody[3];
	const Real aBody[3];
	Shape * shape;
	const double lambda;
	
public:
	OperatorPenalizationMovingFrame(double dt, Real uSolid, Real vSolid, Real wSolid, Real aBody[3], Shape * shape, double lambda) : dt(dt), uBody{uSolid,vSolid,wSolid}, aBody{aBody[0],aBody[1],aBody[2]}, shape(shape), lambda(lambda) {}
	~OperatorPenalizationMovingFrame() {}
	
	void operator()(const BlockInfo& info, FluidBlock& block) const
	{
		Real com[3];
		shape->getCenterOfMass(com);
		
		Real omegaBody[3];
		shape->getAngularVelocity(omegaBody);
		
		// this implementation considers that the Euler updates has already happened
		// do we need a finite state machine coordinating operators?
		for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
			for(int iy=0; iy<FluidBlock::sizeY; ++iy)
				for(int ix=0; ix<FluidBlock::sizeX; ++ix)
				{
					Real p[3], r[3];
					info.pos(p,ix,iy,iz);
					r[0] = p[0] - com[0];
					r[1] = p[1] - com[1];
					r[2] = p[2] - com[2];
					
					Real urBody[3] = {	omegaBody[1]*r[2]-omegaBody[2]*r[1],
						omegaBody[2]*r[0]-omegaBody[0]*r[2],
						omegaBody[0]*r[1]-omegaBody[1]*r[0] };
					
					Real dtChiLambda = dt * lambda * block(ix,iy,iz).chi;
					
					//if (isnan(uBody[0]) || isnan(urBody[0]))
					//	cout << block(ix,iy,iz).u << " " << dt << " " << lambda << " " << block(ix,iy,iz).chi << " " << dtChiLambda << " " << uBody[0] << " " << urBody[0] << endl;
					
					/*
					block(ix,iy,iz).u = (block(ix,iy,iz).u - dt*aBody[0] + dtChiLambda * (uBody[0]+urBody[0])) / (1. + dtChiLambda);
					block(ix,iy,iz).v = (block(ix,iy,iz).v - dt*aBody[1] + dtChiLambda * (uBody[1]+urBody[1])) / (1. + dtChiLambda);
					block(ix,iy,iz).w = (block(ix,iy,iz).w - dt*aBody[2] + dtChiLambda * (uBody[2]+urBody[2])) / (1. + dtChiLambda);
					/*/
					block(ix,iy,iz).u = (block(ix,iy,iz).u + dtChiLambda * (uBody[0]+urBody[0])) / (1. + dtChiLambda);
					block(ix,iy,iz).v = (block(ix,iy,iz).v + dtChiLambda * (uBody[1]+urBody[1])) / (1. + dtChiLambda);
					block(ix,iy,iz).w = (block(ix,iy,iz).w + dtChiLambda * (uBody[2]+urBody[2])) / (1. + dtChiLambda);
					//*/

					
					Real chi = shape->chi(p, info.h_gridpoint);
					block(ix,iy,iz).chi = chi;
					block(ix,iy,iz).rho = shape->rho(p, info.h_gridpoint, chi);
				}
	}
};

#endif
