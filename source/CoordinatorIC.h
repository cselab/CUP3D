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
#include "OperatorIC.h"
#include "Shape.h"

class CoordinatorIC : public GenericCoordinator
{
protected:
	Shape * shape;
	const double uinf;
	
public:
	CoordinatorIC(Shape * shape, const double uinf, FluidGrid * grid) : GenericCoordinator(grid), shape(shape), uinf(uinf)
	{
	}
	
	void operator()(const double dt)
	{
		BlockInfo * ary = &vInfo.front();
		const int N = vInfo.size();
		
		const double dh = vInfo[0].h_gridpoint;
		cout << dh << endl;
		
		Real cx = 0;
		Real cy = 0;
		Real cz = 0;
		Real vol = 0;
		Real gcx = 0;
		Real gcy = 0;
		Real gcz = 0;
		Real gvol = 0;
		Real com[3];
		//shape->getCenterOfMass(com);
		
#pragma omp parallel for reduction(+:cx) reduction(+:cy) reduction(+:cz) reduction(+:vol)
		for(int i=0; i<(int)vInfo.size(); i++)
		{
			BlockInfo info = vInfo[i];
			FluidBlock& b = *(FluidBlock*)info.ptrBlock;
			
			for(int iz=0; iz<FluidBlock::sizeZ; iz++)
				for(int iy=0; iy<FluidBlock::sizeY; iy++)
					for(int ix=0; ix<FluidBlock::sizeX; ix++)
					{
						Real p[3];
						info.pos(p, ix, iy, iz);
						
						b(ix,iy,iz).u = uinf;
						b(ix,iy,iz).v = 0;
						b(ix,iy,iz).w = 0;
						b(ix,iy,iz).chi = shape->chi(p, dh);
						
						// assume fluid with density 1
						b(ix,iy,iz).rho = shape->rho(p, dh, b(ix,iy,iz).chi);
						
						b(ix,iy,iz).p = 0;
						b(ix,iy,iz).divU = 0;
						b(ix,iy,iz).pOld = 0;
						
						b(ix,iy,iz).tmpU = 0;
						b(ix,iy,iz).tmpV = 0;
						b(ix,iy,iz).tmpW = 0;
						b(ix,iy,iz).tmp  = 0;
						
						const Real rhochi = b(ix,iy,iz).rho * b(ix,iy,iz).chi;
						cx += p[0] * rhochi * dh*dh*dh;
						cy += p[1] * rhochi * dh*dh*dh;
						cz += p[2] * rhochi * dh*dh*dh;
						vol += rhochi * dh*dh*dh;
					}
		}
		
		MPI::COMM_WORLD.Allreduce(&cx, &gcx, 1, MPI::DOUBLE, MPI::SUM);
		MPI::COMM_WORLD.Allreduce(&cy, &gcy, 1, MPI::DOUBLE, MPI::SUM);
		MPI::COMM_WORLD.Allreduce(&cz, &gcz, 1, MPI::DOUBLE, MPI::SUM);
		MPI::COMM_WORLD.Allreduce(&vol, &gvol, 1, MPI::DOUBLE, MPI::SUM);
		
		gcx /= gvol;
		gcy /= gvol;
		gcz /= gvol;
		com[0] = gcx;
		com[1] = gcy;
		com[2] = gcz;
		shape->setCenterOfMass(com);
		shape->setCentroid(com);
		
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
	CoordinatorIC_RT(const double rhoS, FluidGrid * grid) : GenericCoordinator(grid), rhoS(rhoS)
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
