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
	CoordinatorBodyVelocities(Real * uBody, Real * vBody, Real * wBody, Real * lambda, Shape * shape, Real * uFlowMax, FluidGridMPI * grid) : GenericCoordinator(grid), uBody(uBody), vBody(vBody), wBody(wBody), lambda(lambda), shape(shape), uFlowMax(uFlowMax)
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
		double J0 = 0;
		double J1 = 0;
		double J2 = 0;
		double J3 = 0;
		double J4 = 0;
		double J5 = 0;
		Real maxU = 0;
		
		double massG = 0;
		double uG = 0;
		double vG = 0;
		double wG = 0;
		double dtdtxG = 0;
		double dtdtyG = 0;
		double dtdtzG = 0;
		double J0G = 0;
		double J1G = 0;
		double J2G = 0;
		double J3G = 0;
		double J4G = 0;
		double J5G = 0;
		Real maxUG = 0;
		
		Real com[3];
		shape->getCenterOfMass(com);
		
#pragma omp parallel for schedule(static) reduction(+:u) reduction(+:v) reduction(+:w), reduction(+:mass)
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
					}
		}
		
		MPI::COMM_WORLD.Allreduce(&u, &uG, 1, MPI::DOUBLE, MPI::SUM);
		MPI::COMM_WORLD.Allreduce(&v, &vG, 1, MPI::DOUBLE, MPI::SUM);
		MPI::COMM_WORLD.Allreduce(&w, &wG, 1, MPI::DOUBLE, MPI::SUM);
		MPI::COMM_WORLD.Allreduce(&mass, &massG, 1, MPI::DOUBLE, MPI::SUM);
		
		*uBody = uG / massG;
		*vBody = vG / massG;
		*wBody = wG / massG;
		
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
						cp[0] = p[1] * (wLocal-*wBody) - p[2] * (vLocal-*vBody); // does this only use the translational component? not from the equations, double check
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
		
		MPI::COMM_WORLD.Allreduce(&dtdtx, &dtdtxG, 1, MPI::DOUBLE, MPI::SUM);
		MPI::COMM_WORLD.Allreduce(&dtdty, &dtdtyG, 1, MPI::DOUBLE, MPI::SUM);
		MPI::COMM_WORLD.Allreduce(&dtdtz, &dtdtzG, 1, MPI::DOUBLE, MPI::SUM);
		MPI::COMM_WORLD.Allreduce(&J0, &J0G, 1, MPI::DOUBLE, MPI::SUM);
		MPI::COMM_WORLD.Allreduce(&J1, &J1G, 1, MPI::DOUBLE, MPI::SUM);
		MPI::COMM_WORLD.Allreduce(&J2, &J2G, 1, MPI::DOUBLE, MPI::SUM);
		MPI::COMM_WORLD.Allreduce(&J3, &J3G, 1, MPI::DOUBLE, MPI::SUM);
		MPI::COMM_WORLD.Allreduce(&J4, &J4G, 1, MPI::DOUBLE, MPI::SUM);
		MPI::COMM_WORLD.Allreduce(&J5, &J5G, 1, MPI::DOUBLE, MPI::SUM);
		MPI::COMM_WORLD.Allreduce(&maxU, &maxUG, 1, MPI::DOUBLE, MPI::MAX);
		
		*uFlowMax = maxUG;
		
		const Real ub[3] = { *uBody, *vBody, *wBody };
		Real dthetadt[3] = { dtdtxG, dtdtyG, dtdtzG };
		const Real weaken = 1.;
		const double J[6] = { J0G*weaken, J1G*weaken, J2G*weaken, J3G*weaken, J4G*weaken, J5G*weaken };
		
		shape->updatePosition(ub, dthetadt, J, massG, dt);
	}
	
	string getName()
	{
		return "BodyVelocities";
	}
};

class CoordinatorBodyVelocitiesForcedRot : public GenericCoordinator
{
protected:
	Real *lambda;
	Real rhoS;
	Shape *shape;
	
public:
	CoordinatorBodyVelocitiesForcedRot(Real * lambda, Shape * shape, FluidGridMPI * grid) : GenericCoordinator(grid), lambda(lambda), shape(shape)
	{
		//cout << "Not supported yet\n";
		//abort();
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
		/*
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
		*/
		
		Real mass = 1;
		
		const Real ub[3] = { 0,0,0 };
		const Real dthetadt[3] = { 0, .1, 0 };
		const double J[6] = { 1,1,1,0,0,0 };//{ J0, J1, J2, J3, J4, J5 };
		
		shape->updatePosition(ub, dthetadt, J, mass, dt);
	}
	
	string getName()
	{
		return "BodyVelocities";
	}
};


#endif
