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
	Real *uFlowMax;
	Real *lambda;
	Real rhoS;
	Shape *shape;
	
public:
	CoordinatorBodyVelocities(Real * uBody, Real * vBody, Real * wBody, Real * lambda, Shape * shape, Real * uFlowMax, FluidGrid * grid) : GenericCoordinator(grid), uBody(uBody), vBody(vBody), wBody(wBody), lambda(lambda), shape(shape), uFlowMax(uFlowMax)
	{
	}
	
	void operator()(const double dt)
	{
		double mass = 0;
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
		Real maxU = 0;
		Real cx = 0;
		Real cy = 0;
		Real cz = 0;
		Real vol = 0;
		
		Real com[3];
		shape->getCenterOfMass(com);
		
#pragma omp parallel for schedule(static) reduction(+:u) reduction(+:v) reduction(+:w), reduction(+:mass) reduction(+:cx) reduction(+:cy) reduction(+:cz) reduction(+:vol)
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
						Real p[3];
						info.pos(p, ix, iy, iz);
						
						Real uLocal = b(ix,iy,iz).u;
						Real vLocal = b(ix,iy,iz).v;
						Real wLocal = b(ix,iy,iz).w;
						
						const Real rhochi = b(ix,iy,iz).rho * b(ix,iy,iz).chi;
						u += uLocal * rhochi * h3;
						v += vLocal * rhochi * h3;
						w += wLocal * rhochi * h3;
						
						mass += rhochi * h3;
						
						cx += p[0] * rhochi * h3;
						cy += p[1] * rhochi * h3;
						cz += p[2] * rhochi * h3;
						vol += rhochi * h3;
					}
		}
		
		cx /= vol;
		cy /= vol;
		cz /= vol;
		//cout << "\tGeometricCoM:\t" << com[0] << " " << com[1] << " " << com[2] << endl;
		//cout << "\tGridCM:\t" << cx << " " << cy << " " << cz << endl;
		
		*uBody = u / mass;
		*vBody = v / mass;
		*wBody = w / mass;
		
#pragma omp parallel for schedule(static) reduction(+:dtdtx) reduction(+:dtdty) reduction(+:dtdtz) reduction(+:J0) reduction(+:J1) reduction(+:J2) reduction(+:J3) reduction(+:J4) reduction(+:J5) reduction(max:maxU)
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
						Real uLocal = b(ix,iy,iz).u;
						Real vLocal = b(ix,iy,iz).v;
						Real wLocal = b(ix,iy,iz).w;
						
#ifndef _MOVING_FRAME_
						maxU = max(maxU,(Real)abs(uLocal));
						maxU = max(maxU,(Real)abs(vLocal));
						maxU = max(maxU,(Real)abs(wLocal));
#else
						maxU = max(maxU,(Real)abs(uLocal-*uBody));
						maxU = max(maxU,(Real)abs(vLocal-*vBody));
						maxU = max(maxU,(Real)abs(wLocal-*wBody));
#endif
						
						const Real rhochi = b(ix,iy,iz).rho * b(ix,iy,iz).chi;
						Real p[3];
						info.pos(p, ix, iy, iz);
						p[0] -= com[0];
						p[1] -= com[1];
						p[2] -= com[2];
						
						Real cp[3];
						cp[0] = p[1] * (wLocal-*wBody) - p[2] * (vLocal-*vBody);
						cp[1] = p[2] * (uLocal-*uBody) - p[0] * (wLocal-*wBody);
						cp[2] = p[0] * (vLocal-*vBody) - p[1] * (uLocal-*uBody);
						//cp[0] = p[1] * (wLocal) - p[2] * (vLocal);
						//cp[1] = p[2] * (uLocal) - p[0] * (wLocal);
						//cp[2] = p[0] * (vLocal) - p[1] * (uLocal);
						
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
		
		*uFlowMax = maxU;
		
		const Real ub[3] = { *uBody, *vBody, *wBody };
		Real dthetadt[3] = { dtdtx, dtdty, dtdtz };
		const Real weaken = .001;
		const Real J[6] = { J0*weaken, J1*weaken, J2*weaken, J3*weaken, J4*weaken, J5*weaken };
		
		shape->updatePosition(ub, dthetadt, J, mass, dt);
	}
	
	string getName()
	{
		return "BodyVelocities";
	}
};

class CoordinatorBodyVelocitiesForcedRot : public GenericCoordinator
{
protected:
	Real *omegaBodyX, *omegaBodyY, *omegaBodyZ;
	Real *lambda;
	Real rhoS;
	Shape *shape;
	
public:
	CoordinatorBodyVelocitiesForcedRot(Real * omegaBodyX, Real * omegaBodyY, Real * omegaBodyZ, Real * lambda, Shape * shape, FluidGrid * grid) : GenericCoordinator(grid), omegaBodyX(omegaBodyX), omegaBodyY(omegaBodyY), omegaBodyZ(omegaBodyZ), lambda(lambda), shape(shape)
	{
	}
	
	void operator()(const double dt)
	{
		Real maxU = 0;
		const int N = vInfo.size();
		Real J0 = 0;
		Real J1 = 0;
		Real J2 = 0;
		Real J3 = 0;
		Real J4 = 0;
		Real J5 = 0;
		Real dtdtx = 0;
		Real dtdty = 0;
		Real dtdtz = 0;
		
		
		Real com[3];
		shape->getCenterOfMass(com);
		
#pragma omp parallel for schedule(static) reduction(+:J0) reduction(+:J1) reduction(+:J2) reduction(+:J3) reduction(+:J4) reduction(+:J5) reduction(+:dtdtx) reduction(+:dtdty) reduction(+:dtdtz)
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
						Real p[3];
						info.pos(p, ix, iy, iz);
						
						const Real rhochi = b(ix,iy,iz).rho * b(ix,iy,iz).chi;
						p[0] -= com[0];
						p[1] -= com[1];
						p[2] -= com[2];
						
						const Real uLocal = -p[2] * .2 * M_PI;
						const Real vLocal = 0;
						const Real wLocal =  p[0] * .2 * M_PI;
						
						Real cp[3];
						cp[0] = p[1]* wLocal - p[2]* vLocal;
						cp[1] = p[2]* uLocal - p[0]* wLocal;
						cp[2] = p[0]* vLocal - p[1]* uLocal;
						
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
		
		Real mass = 1;
		
		const Real ub[3] = { 0,0,0 };
		const Real dthetadt[3] = { 0, 0, .1 };
		const Real J[6] = { J0, J1, J2, J3, J4, J5 };
		
		shape->updatePosition(ub, dthetadt, J, mass, dt);
	}
	
	string getName()
	{
		return "BodyVelocities";
	}
};


#endif
