//
//  CoordinatorIC.h
//  CubismUP_3D
//
//  Created by Christian Conti on 3/27/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef CubismUP_3D_CoordinatorIC_h
#define CubismUP_3D_CoordinatorIC_h

#include "GenericCoordinator.h"
#include "GenericOperator.h"

class OperatorIC : public GenericOperator
{
private:
	const double uinf;

public:
	OperatorIC(const double uinf) : uinf(uinf) {}

	~OperatorIC() {}

	void operator()(const BlockInfo& info, FluidBlock& block) const
	{
		for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
			for(int iy=0; iy<FluidBlock::sizeY; ++iy)
				for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
					block(ix,iy,iz).u = uinf;
					block(ix,iy,iz).v = 0;
					block(ix,iy,iz).w = 0;
					block(ix,iy,iz).chi = 0;

					block(ix,iy,iz).p = 0;
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
public:
	OperatorIC_RT(const double rhoS) {}

	~OperatorIC_RT() {}

	void operator()(const BlockInfo& info, FluidBlock& block) const
	{
		Real d = .25;
		Real a = .05;//.25

		for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
			for(int iy=0; iy<FluidBlock::sizeY; ++iy)
				for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
					Real p[3];
					info.pos(p, ix, iy, iz);

					Real x = p[0] - d*.5;
					Real y = p[1] - d*2.5;
					Real z = p[2] - d*.5;
					Real eta = -.1*a*cos(2*M_PI*x/d);

					//block(ix,iy).rho = 2. + tanh((y-eta)/(.01*d));

					Real eta = -a*.25*cos(2*M_PI*x/d)*cos(2*M_PI*z/d);

					block(ix,iy,iz).u = 0;
					block(ix,iy,iz).v = 0;
					block(ix,iy,iz).w = 0;
					block(ix,iy,iz).chi = .5+.5*tanh((y-eta)/info.h_gridpoint);

					block(ix,iy,iz).p = 0;
					block(ix,iy,iz).pOld = 0;

					block(ix,iy,iz).tmpU = 0;
					block(ix,iy,iz).tmpV = 0;
					block(ix,iy,iz).tmpW = 0;
					block(ix,iy,iz).tmp  = 0;
				}
	}
};

class CoordinatorIC : public GenericCoordinator
{
protected:
	const double uinf;
	
public:
	CoordinatorIC(const double uinf, FluidGridMPI * grid) : GenericCoordinator(grid), uinf(uinf)
	{
	}
	
	void operator()(const double dt)
	{
		BlockInfo * ary = &vInfo.front();
		const int N = vInfo.size();
		
#pragma omp parallel
		{
			OperatorIC kernel(rhoS);
			
#pragma omp for schedule(static)
			for (int i=0; i<N; i++)
				kernel(ary[i], *(FluidBlock*)ary[i].ptrBlock);
		}
		check("IC - end");
	}
	
	string getName()
	{
		return "IC";
	}
};

class CoordinatorIC_RT : public GenericCoordinator
{
protected:
	const double rhoS;
	
public:
	CoordinatorIC_RT(const double rhoS, FluidGridMPI * grid) : GenericCoordinator(grid), rhoS(rhoS)
	{
	}
	
	void operator()(const double dt)
	{
		BlockInfo * ary = &vInfo.front();
		const int N = vInfo.size();
		
#pragma omp parallel
		{
			OperatorIC_RT kernel(rhoS);
			
#pragma omp for schedule(static)
			for (int i=0; i<N; i++)
				kernel(ary[i], *(FluidBlock*)ary[i].ptrBlock);
		}
		
		check("IC - end");
	}
	
	string getName()
	{
		return "IC_RT";
	}
};

#endif
