//
//  CoordinatorBodyVelocities.h
//  CubismUP_3D
//
//  Created by Christian Conti on 3/30/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef CubismUP_3D_CoordinatorComputeBodyVelocities_h
#define CubismUP_3D_CoordinatorComputeBodyVelocities_h

class CoordinatorBodyVelocities : public GenericCoordinator
{
protected:
	Real *uBody, *vBody, *omegaBody;
	Real *lambda;
	Real rhoS;
	
public:
	CoordinatorBodyVelocities(Real * uBody, Real * vBody, Real * omegaBody, Real * lambda, Real rhoS, FluidGrid * grid) : GenericCoordinator(grid), uBody(uBody), vBody(vBody), omegaBody(omegaBody), lambda(lambda), rhoS(rhoS)
	{
	}
	
	void operator()(const double dt)
	{
		double centerTmpX = 0;
		double centerTmpY = 0;
		double mass = 0;
		double volume = 0;
		double u = 0;
		double v = 0;
		double momOfInertia = 0;
		double angularMomentum = 0;
		const int N = vInfo.size();
		
#pragma omp parallel for schedule(static) reduction(+:centerTmpX) reduction(+:centerTmpY) reduction(+:mass)
		for(int i=0; i<N; i++)
		{
			BlockInfo info = vInfo[i];
			FluidBlock& b = *(FluidBlock*)info.ptrBlock;
			
			double h = info.h_gridpoint;
			
			for(int iy=0; iy<FluidBlock::sizeY; ++iy)
				for(int ix=0; ix<FluidBlock::sizeX; ++ix)
				{
					double p[2] = {0,0};
					info.pos(p, ix, iy);
					double rhochi = b(ix,iy).rho * b(ix,iy).chi;
					centerTmpX += p[0] * rhochi;
					centerTmpY += p[1] * rhochi;
					mass += rhochi;
					volume += b(ix,iy).chi;
				}
		}
		
        // needs to be fixed for periodicity
		centerTmpX /= mass;
		centerTmpY /= mass;
		
		//*
#pragma omp parallel for schedule(static) reduction(+:u) reduction(+:v) reduction(+:momOfInertia) reduction(+:angularMomentum)
		for(int i=0; i<N; i++)
		{
			BlockInfo info = vInfo[i];
			FluidBlock& b = *(FluidBlock*)info.ptrBlock;
			
			double h = info.h_gridpoint;
			
			for(int iy=0; iy<FluidBlock::sizeY; ++iy)
				for(int ix=0; ix<FluidBlock::sizeX; ++ix)
				{
					double p[2] = {0,0};
					info.pos(p, ix, iy);
#ifndef _AVGU_
					double rhochi = b(ix,iy).rho * b(ix,iy).chi;
					u += b(ix,iy).u * rhochi;
					v += b(ix,iy).v * rhochi;
					momOfInertia    += rhochi * ((p[0]-centerTmpX)*(p[0]-centerTmpX) + (p[1]-centerTmpY)*(p[1]-centerTmpY));
					angularMomentum += rhochi * ((p[0]-centerTmpX)*b(ix,iy).v        - (p[1]-centerTmpY)*b(ix,iy).u);
#else
					double chi = b(ix,iy).chi;
					u += b(ix,iy).u * chi;
					v += b(ix,iy).v * chi;
					momOfInertia    += chi * ((p[0]-centerTmpX)*(p[0]-centerTmpX) + (p[1]-centerTmpY)*(p[1]-centerTmpY));
					angularMomentum += chi * ((p[0]-centerTmpX)*b(ix,iy).v        - (p[1]-centerTmpY)*b(ix,iy).u);
#endif
				}
		}
	
#ifndef _AVGU_
		*uBody = u / mass;
		*vBody = v / mass;
#else
		*uBody = u / volume;
		*vBody = v / volume;
#endif
		*omegaBody = angularMomentum / momOfInertia;
		
		/*/
		u=0;
		v=0;
#pragma omp parallel for schedule(static) reduction(+:u) reduction(+:v)
		for(int i=0; i<N; i++)
		{
		 BlockInfo info = vInfo[i];
		 FluidBlock& b = *(FluidBlock*)info.ptrBlock;
		 
		 Real h = info.h_gridpoint;
		 
		 for(int iy=0; iy<FluidBlock::sizeY; ++iy)
			for(int ix=0; ix<FluidBlock::sizeX; ++ix)
			{
				u += (b(ix,iy).u-(*uBody)) * b(ix,iy).chi;
				v += (b(ix,iy).v-(*vBody)) * b(ix,iy).chi;
			}
		}
		
		*uBody += dt*u*(*lambda) / mass;
		*vBody += dt*v*(*lambda) / mass;
		
		cout << "vBody is " << *vBody << endl;
		//*/
	}
	
	string getName()
	{
		return "BodyVelocities";
	}
};


#endif
