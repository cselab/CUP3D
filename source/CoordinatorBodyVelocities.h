//
//  CoordinatorBodyVelocities.h
//  CubismUP_3D
//
//  Created by Christian Conti on 3/30/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef CubismUP_3D_CoordinatorComputeBodyVelocities_h
#define CubismUP_3D_CoordinatorComputeBodyVelocities_h

#include "GenericCoordinator.h"
#include "Shape.h"

class CoordinatorBodyVelocities : public GenericCoordinator
{
protected:
	Real *uBody, *vBody, *wBody;
	Real *lambda;
	Real rhoS;
	Shape *shape;
	
public:
	CoordinatorBodyVelocities(Real * uBody, Real * vBody, Real * wBody, Real * lambda, Shape * shape, FluidGrid * grid) : GenericCoordinator(grid), uBody(uBody), vBody(vBody), wBody(wBody), lambda(lambda), shape(shape)
	{
	}
	
	void operator()(const double dt)
	{
		double centerTmpX = 0;
		double centerTmpY = 0;
		double centerTmpZ = 0;
		double mass = 0;
		double volume = 0;
		double u = 0;
		double v = 0;
		double w = 0;
		double dtdtx = 0;
		double dtdty = 0;
		double dtdtz = 0;
		const int N = vInfo.size();
		Real J0 = 0;
		Real J1 = 0;
		Real J2 = 0;
		Real J3 = 0;
		Real J4 = 0;
		Real J5 = 0;
		
#pragma omp parallel for schedule(static) reduction(+:centerTmpX) reduction(+:centerTmpY) reduction(+:centerTmpZ) reduction(+:mass)
		for(int i=0; i<N; i++)
		{
			BlockInfo info = vInfo[i];
			FluidBlock& b = *(FluidBlock*)info.ptrBlock;
			
			const Real h = info.h_gridpoint;
			const Real h3 = h*h*h;
			
			for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
				for(int iy=0; iy<FluidBlock::sizeY; ++iy)
					for(int ix=0; ix<FluidBlock::sizeX; ++ix)
					{
						double p[3] = {0,0,0};
						info.pos(p, ix, iy, iz);
						double rhochi = b(ix,iy,iz).rho * b(ix,iy,iz).chi;
						centerTmpX += p[0] * rhochi * h3;
						centerTmpY += p[1] * rhochi * h3;
						centerTmpZ += p[2] * rhochi * h3;
						mass += rhochi * h3;
						volume += b(ix,iy,iz).chi * h3;
					}
		}
		
		// needs to be fixed for periodicity
		centerTmpX /= mass;
		centerTmpY /= mass;
		centerTmpZ /= mass;
		
		//Real com[3] = {centerTmpX,centerTmpY,centerTmpZ};
		//shape->setCenterOfMass(com);
		
		//*
#pragma omp parallel for schedule(static) reduction(+:u) reduction(+:v) reduction(+:w) reduction(+:dtdtx) reduction(+:dtdty) reduction(+:dtdtz) reduction(+:J0) reduction(+:J1) reduction(+:J2) reduction(+:J3) reduction(+:J4) reduction(+:J5)
		for(int i=0; i<N; i++)
		{
			BlockInfo info = vInfo[i];
			FluidBlock& b = *(FluidBlock*)info.ptrBlock;
			
			const Real h = info.h_gridpoint;
			const Real h3 = h*h*h;
			
			for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
				for(int iy=0; iy<FluidBlock::sizeY; ++iy)
					for(int ix=0; ix<FluidBlock::sizeX; ++ix)
					{
#ifndef _AVGU_
						const Real rhochi = b(ix,iy,iz).rho * b(ix,iy,iz).chi;
						u += b(ix,iy,iz).u * rhochi * h3;
						v += b(ix,iy,iz).v * rhochi * h3;
						w += b(ix,iy,iz).w * rhochi * h3;
#else
						const Real chi = b(ix,iy,iz).chi;
						u += b(ix,iy,iz).u * chi * h3;
						v += b(ix,iy,iz).v * chi * h3;
						w += b(ix,iy,iz).w * chi * h3;
#endif
						double p[3] = {0,0,0};
						info.pos(p, ix, iy, iz);
						p[0] -= centerTmpX;
						p[1] -= centerTmpY;
						p[2] -= centerTmpZ;
						
						Real cp[3];
						cp[0] = p[1]* b(ix,iy,iz).w - p[2]* b(ix,iy,iz).v;
						cp[1] = p[2]* b(ix,iy,iz).u - p[0]* b(ix,iy,iz).w;
						cp[2] = p[0]* b(ix,iy,iz).v - p[1]* b(ix,iy,iz).u;
						
						dtdtx += cp[0] * rhochi * h3;
						dtdty += cp[1] * rhochi * h3;
						dtdtz += cp[2] * rhochi * h3;
						
						J0 += rhochi * (p[1]*p[1] + p[2]*p[2]) * h3; //       y^2 + z^2
						J1 += rhochi * (p[0]*p[0] + p[2]*p[2]) * h3; // x^2 +     + z^2
						J2 += rhochi * (p[0]*p[0] + p[1]*p[1]) * h3; // x^2 + y^2
						J3 -= rhochi * p[0] * p[1] * h3; // xy
						J4 -= rhochi * p[0] * p[2] * h3; // xz
						J5 -= rhochi * p[1] * p[2] * h3; // yz
					}
		}
		
#ifndef _AVGU_
		*uBody = u / mass;
		*vBody = v / mass;
		*wBody = w / mass;
#else
		*uBody = u / volume;
		*vBody = v / volume;
		*wBody = w / volume;
#endif
		
		const Real ub[3] = { *uBody, *vBody, *wBody };
		const Real dthetadt[3] = { dtdtx, dtdty, dtdtz };
		const Real J[6] = { J0, J1, J2, J3, J4, J5 };
		//const Real dthetadt[3] = { 0, 0, 2*M_PI };
		//const Real J[6] = { 1, 1, 1, 0, 0, 0 };
		
		//cout << "\tMass:\t" << mass << endl;
		//cout << "\tCoM:\t" << centerTmpX << " " << centerTmpY << " " << centerTmpZ << endl;
		//cout << "\tL:\t" << dtdtx << " " << dtdty << " " << dtdtz << endl;
		//cout << "\tJ:\t" << J0 << " " << J1 << " " << J2 << " " << J3 << " " << J4 << " " << J5 << endl;
		shape->updatePosition(ub, dthetadt, J, mass, dt);
	}
	
	string getName()
	{
		return "BodyVelocities";
	}
};


#endif
