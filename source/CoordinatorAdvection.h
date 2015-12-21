//
//  CoordinatorAdvection.h
//  CubismUP_3D
//
//  Created by Christian Conti on 3/27/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef CubismUP_3D_CoordinatorAdvection_h
#define CubismUP_3D_CoordinatorAdvection_h

#include "GenericCoordinator.h"
#include "OperatorAdvection.h"
#include <cmath>

template <typename Lab>
class CoordinatorAdvection : public GenericCoordinator
{
protected:
	Real *uBody, *vBody, *wBody;
#ifdef _MULTIPHASE_
	Real rhoS;
#endif
	
	inline void update()
	{
		const int N = vInfo.size();
		
#pragma omp parallel for schedule(static)
		for(int i=0; i<N; i++)
		{
			BlockInfo info = vInfo[i];
			FluidBlock& b = *(FluidBlock*)info.ptrBlock;
			
			for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
				for(int iy=0; iy<FluidBlock::sizeY; ++iy)
					for(int ix=0; ix<FluidBlock::sizeX; ++ix)
					{
						b(ix,iy,iz).u = b(ix,iy,iz).tmpU;
						b(ix,iy,iz).v = b(ix,iy,iz).tmpV;
						b(ix,iy,iz).w = b(ix,iy,iz).tmpW;
#ifdef _MULTIPHASE_
						//b(ix,iy).chi = b(ix,iy).tmp;
						//b(ix,iy).rho = b(ix,iy).chi * rhoS + (1-b(ix,iy).chi);
						
						// threshold density
						//b(ix,iy).rho = b(ix,iy).tmp;
						Real density = min(max(b(ix,iy,iz).tmp,min((Real)1.,rhoS)),max((Real)1.,rhoS));
						b(ix,iy,iz).rho = density;
#endif // _MULTIPHASE_
					}
		}
	}
	
	inline void advect(const double dt, const int stage)
	{
		OperatorAdvectionUpwind3rdOrder kernel(dt,uBody,vBody,wBody,stage);
		compute(kernel);
		/*
		BlockInfo * ary = &vInfo.front();
		const int N = vInfo.size();
		
#pragma omp parallel
		{
			OperatorAdvectionUpwind3rdOrder kernel(dt,uBody,vBody,wBody,stage);
			
			Lab lab;
			lab.prepare(*grid, kernel.stencil_start, kernel.stencil_end, false);
			
#pragma omp for schedule(static)
			for (int i=0; i<N; i++)
			{
				lab.load(ary[i], 0);
				kernel(lab, ary[i], *(FluidBlock*)ary[i].ptrBlock);
			}
		}
		 */
	}
	
public:
#ifndef _MULTIPHASE_
	CoordinatorAdvection(Real * uBody, Real * vBody, Real * wBody, FluidGrid * grid) : GenericCoordinator(grid), uBody(uBody), vBody(vBody), wBody(wBody)
#else
	CoordinatorAdvection(Real * uBody, Real * vBody, Real * wBody, FluidGrid * grid, Real rhoS) : GenericCoordinator(grid), uBody(uBody), vBody(vBody), wBody(wBody), rhoS(rhoS)
#endif
	{
	}
	
#ifndef _MULTIPHASE_
	CoordinatorAdvection(FluidGrid * grid) : GenericCoordinator(grid), uBody(NULL), vBody(NULL)
#else
	CoordinatorAdvection(FluidGrid * grid, Real rhoS) : GenericCoordinator(grid), uBody(NULL), vBody(NULL), wBody(NULL), rhoS(rhoS)
#endif
	{
	}
	
	void operator()(const double dt)
	{
		check("advection - start");
		
		advect(dt,0);
#ifdef _RK2_
		advect(dt,1);
#endif
		update();
		
		check("advection - end");
	}
	
	string getName()
	{
		return "Advection";
	}
};

template <typename Lab>
class CoordinatorTransport : public GenericCoordinator
{
protected:
	inline void reset()
	{
		const int N = vInfo.size();
		
#pragma omp parallel for schedule(static)
		for(int i=0; i<N; i++)
		{
			BlockInfo info = vInfo[i];
			FluidBlock& b = *(FluidBlock*)info.ptrBlock;
			
			for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
				for(int iy=0; iy<FluidBlock::sizeY; ++iy)
					for(int ix=0; ix<FluidBlock::sizeX; ++ix)
						b(ix,iy,iz).tmp = 0;
		}
	};
	
	inline void update()
	{
		const int N = vInfo.size();
		
#pragma omp parallel for schedule(static)
		for(int i=0; i<N; i++)
		{
			BlockInfo info = vInfo[i];
			FluidBlock& b = *(FluidBlock*)info.ptrBlock;
			
			for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
				for(int iy=0; iy<FluidBlock::sizeY; ++iy)
					for(int ix=0; ix<FluidBlock::sizeX; ++ix)
						b(ix,iy,iz).rho = b(ix,iy,iz).tmp;
		}
	}
	
public:
	CoordinatorTransport(FluidGrid * grid) : GenericCoordinator(grid)
	{
	}
	
	void operator()(const double dt)
	{
		BlockInfo * ary = &vInfo.front();
		const int N = vInfo.size();
		
		reset();
		
#pragma omp parallel
		{
			OperatorTransportUpwind3rdOrder kernel(dt,0);
			
			Lab lab;
			lab.prepare(*grid, kernel.stencil_start, kernel.stencil_end, false);
			
#pragma omp for schedule(static)
			for (int i=0; i<N; i++)
			{
				lab.load(ary[i], 0);
				kernel(lab, ary[i], *(FluidBlock*)ary[i].ptrBlock);
			}
		}
		
#ifdef _RK2_
#pragma omp parallel
		{
			OperatorTransportUpwind3rdOrder kernel(dt,1);
			
			Lab lab;
			lab.prepare(*grid, kernel.stencil_start, kernel.stencil_end, true);
			
#pragma omp for schedule(static)
			for (int i=0; i<N; i++)
			{
				lab.load(ary[i], 0);
				kernel(lab, ary[i], *(FluidBlock*)ary[i].ptrBlock);
			}
		}
#endif // _RK2_
		
		update();
	}
	
	string getName()
	{
		return "Transport";
	}
};

template <typename Lab>
class CoordinatorTransportTimeTest : public GenericCoordinator
{
protected:
	double time;
	
	double _analyticalRHS(double px, double py, double pz, double t)
	{
		return 8 * M_PI * cos((px+t) * 8. * M_PI);
	}
	
	inline void reset()
	{
		const int N = vInfo.size();
		
#pragma omp parallel for schedule(static)
		for(int i=0; i<N; i++)
		{
			BlockInfo info = vInfo[i];
			FluidBlock& b = *(FluidBlock*)info.ptrBlock;
			
			for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
				for(int iy=0; iy<FluidBlock::sizeY; ++iy)
					for(int ix=0; ix<FluidBlock::sizeX; ++ix)
						b(ix,iy,iz).tmp = 0;
		}
	};
	
	inline void update()
	{
		const int N = vInfo.size();
		
#pragma omp parallel for schedule(static)
		for(int i=0; i<N; i++)
		{
			BlockInfo info = vInfo[i];
			FluidBlock& b = *(FluidBlock*)info.ptrBlock;
			
			for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
				for(int iy=0; iy<FluidBlock::sizeY; ++iy)
					for(int ix=0; ix<FluidBlock::sizeX; ++ix)
						b(ix,iy,iz).rho = b(ix,iy,iz).tmp;
		}
	}
	
	inline void advect(const double dt, const int stage)
	{
		BlockInfo * ary = &vInfo.front();
		const int N = vInfo.size();
#pragma omp parallel
		{
			OperatorTransportUpwind3rdOrder kernel(dt,stage);
			
			Lab lab;
			lab.prepare(*grid, kernel.stencil_start, kernel.stencil_end, false);
			
#pragma omp for schedule(static)
			for (int i=0; i<N; i++)
			{
				lab.load(ary[i], 0);
				kernel(lab, ary[i], *(FluidBlock*)ary[i].ptrBlock);
			}
		}
	}
	
public:
	CoordinatorTransportTimeTest(FluidGrid * grid) : GenericCoordinator(grid)
	{
	}
	
	void operator()(const double dt)
	{
		check("advection - start");
		
		reset();
		advect(dt,0);
#ifdef _RK2_
		advect(dt,1);
#endif
		update();
		time+=dt;
		
		check("advection - end");
	}
	
	string getName()
	{
		return "Transport";
	}
};

#endif
