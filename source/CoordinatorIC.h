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
	CoordinatorIC(Shape * shape, const double uinf, FluidGridMPI * grid) : GenericCoordinator(grid), shape(shape), uinf(uinf)
	{
	}
	
	void operator()(const double dt)
	{
		BlockInfo * ary = &vInfo.front();
		const int N = vInfo.size();
		
		const double dh = vInfo[0].h_gridpoint;
		
		double cx = 0;
		double cy = 0;
		double cz = 0;
		double vol = 0;
		double gcx = 0;
		double gcy = 0;
		double gcz = 0;
		double gvol = 0;
		Real com[3];
		//shape->getCenterOfMass(com);
		
		double J0 = 0;
		double J1 = 0;
		double J2 = 0;
		double J3 = 0;
		double J4 = 0;
		double J5 = 0;
		double J0G = 0;
		double J1G = 0;
		double J2G = 0;
		double J3G = 0;
		double J4G = 0;
		double J5G = 0;
		
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
						cx += p[0] * rhochi;
						cy += p[1] * rhochi;
						cz += p[2] * rhochi;
						vol += rhochi;
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
		com[1] = gcy*_SY_;//*1.02;//*1.05;
		com[2] = gcz*_SZ_;//*0.98;//*0.95;
		shape->setCenterOfMass(com);
		shape->setCentroid(com);
		shape->setMass(gvol);
		
		cout << "Center of mass (after IC) set to " << com[0] << " " << com[1] << " " << com[2] << endl;
		
		
#pragma omp parallel for reduction(+:J0) reduction(+:J1) reduction(+:J2) reduction(+:J3) reduction(+:J4) reduction(+:J5)
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
						p[0] -= com[0];
						p[1] -= com[1];
						p[2] -= com[2];
						
						const Real rhochi = b(ix,iy,iz).rho * b(ix,iy,iz).chi;
						
						J0 += rhochi * (p[1]*p[1] + p[2]*p[2]); //       y^2 + z^2
						J1 += rhochi * (p[0]*p[0] + p[2]*p[2]); // x^2 +     + z^2
						J2 += rhochi * (p[0]*p[0] + p[1]*p[1]); // x^2 + y^2
						J3 -= rhochi * p[0] * p[1]; // xy
						J4 -= rhochi * p[0] * p[2]; // xz
						J5 -= rhochi * p[1] * p[2]; // yz
					}
		}
		
		MPI::COMM_WORLD.Allreduce(&J0, &J0G, 1, MPI::DOUBLE, MPI::SUM);
		MPI::COMM_WORLD.Allreduce(&J1, &J1G, 1, MPI::DOUBLE, MPI::SUM);
		MPI::COMM_WORLD.Allreduce(&J2, &J2G, 1, MPI::DOUBLE, MPI::SUM);
		MPI::COMM_WORLD.Allreduce(&J3, &J3G, 1, MPI::DOUBLE, MPI::SUM);
		MPI::COMM_WORLD.Allreduce(&J4, &J4G, 1, MPI::DOUBLE, MPI::SUM);
		MPI::COMM_WORLD.Allreduce(&J5, &J5G, 1, MPI::DOUBLE, MPI::SUM);
		
		
		//const double h = vInfo[0].h_gridpoint;
		//const double h3 = h*h*h;
		
		const double J[6] = { J0G, J1G, J2G, J3G, J4G, J5G };
		
		shape->setInertiaMatrix0(J);
		
		cout << "Inertia matrix (after IC, without h3) set to " << J[0] << " " << J[1] << " " << J[2] << " " << J[3] << " " << J[4] << " " << J[5] << endl;
		
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
