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
	const Real uinf;

public:
	OperatorIC(const Real uinf) : uinf(uinf) {}

	~OperatorIC() {}

	void operator()(const BlockInfo& info, FluidBlock& block) const
	{
		for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
			for(int iy=0; iy<FluidBlock::sizeY; ++iy)
				for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
					/*
					Real p[3];
					info.pos(p, ix, iy, iz);
					const Real x = 2*M_PI*p[0];
					const Real y = 4*M_PI*p[1];
					const Real z = 4*M_PI*p[2];
					 */
					block(ix,iy,iz).u = 0;//.05*cos(x+0.00*M_PI)*cos(y+0.75*M_PI)*cos(z+1.50*M_PI);
					block(ix,iy,iz).v = 0;//.05*cos(x+0.25*M_PI)*cos(y+1.00*M_PI)*cos(z+1.75*M_PI);
					block(ix,iy,iz).w = 0;//.05*cos(x+0.50*M_PI)*cos(y+1.25*M_PI)*cos(z+2.00*M_PI);
					block(ix,iy,iz).chi = 0;
					block(ix,iy,iz).p = 0;
					block(ix,iy,iz).tmpU = 0;
					block(ix,iy,iz).tmpV = 0;
					block(ix,iy,iz).tmpW = 0;
				}
	}
};

class OperatorIC_RT : public GenericOperator
{
public:
	OperatorIC_RT(const Real rhoS) {}

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
					//Real eta = -.1*a*cos(2*M_PI*x/d);
					Real eta = -a*.25*cos(2*M_PI*x/d)*cos(2*M_PI*z/d);

					block(ix,iy,iz).u = 0;
					block(ix,iy,iz).v = 0;
					block(ix,iy,iz).w = 0;
					block(ix,iy,iz).chi = .5+.5*tanh((y-eta)/info.h_gridpoint);

					block(ix,iy,iz).p = 0;

					block(ix,iy,iz).tmpU = 0;
					block(ix,iy,iz).tmpV = 0;
					block(ix,iy,iz).tmpW = 0;
				}
	}
};

class CoordinatorIC : public GenericCoordinator
{
public:
	CoordinatorIC(FluidGridMPI * grid) : GenericCoordinator(grid)
	{
	}
	
	void operator()(const Real dt)
	{
		const int N = vInfo.size();
		
#pragma omp parallel
		{
			OperatorIC kernel(0.);
#pragma omp for schedule(static)
			for (int i=0; i<N; i++) {
				BlockInfo info = vInfo[i];
				FluidBlock& b = *(FluidBlock*)info.ptrBlock;
				kernel(info, b);
			}
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
public:
	CoordinatorIC_RT(FluidGridMPI * grid) : GenericCoordinator(grid)
	{
	}
	
	void operator()(const Real dt)
	{
		const int N = vInfo.size();
		
#pragma omp parallel
		{
			OperatorIC_RT kernel(0.);
#pragma omp for schedule(static)
			for (int i=0; i<N; i++) {
				BlockInfo info = vInfo[i];
				FluidBlock& b = *(FluidBlock*)info.ptrBlock;
				kernel(info, b);
			}
		}
		
		check("IC - end");
	}
	
	string getName()
	{
		return "IC_RT";
	}
};

#endif
