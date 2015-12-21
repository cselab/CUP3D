//
//  CoordinatorDiffusion.h
//  CubismUP_3D
//
//  Created by Christian Conti on 3/27/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef CubismUP_3D_CoordinatorDiffusion_h
#define CubismUP_3D_CoordinatorDiffusion_h

#include "GenericCoordinator.h"
#include "OperatorDiffusion.h"

template <typename Lab>
class CoordinatorDiffusion : public GenericCoordinator
{
protected:
	const double coeff;
	Real *uBody, *vBody, *wBody;
	
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
#ifdef _DENSITYDIFF_
						b(ix,iy,iz).rho = b(ix,iy,iz).tmp;
#endif
					}
		}
	}
	
	inline void diffuse(const double dt, const int stage)
	{
		OperatorDiffusion kernel(dt, coeff, stage);
		compute(kernel);
	}
	
public:
	CoordinatorDiffusion(const double coeff, Real * uBody, Real * vBody, Real * wBody, FluidGrid * grid) : GenericCoordinator(grid), coeff(coeff), uBody(uBody), vBody(vBody), wBody(wBody)
	{
	}
	
	CoordinatorDiffusion(const double coeff, FluidGrid * grid) : GenericCoordinator(grid), coeff(coeff), uBody(NULL), vBody(NULL), wBody(NULL)
	{
	}
	
	void operator()(const double dt)
	{
		check("diffusion - start");
		
		diffuse(dt,0);
#ifdef _RK2_
		diffuse(dt,1);
#endif
		update();
		
		check("diffusion - end");
	}
	
	string getName()
	{
		return "Diffusion";
	}
};


class CoordinatorDiffusionTimeTest : public GenericCoordinator
{
protected:
	const double coeff;
	double time;
	const double freq;
	
	double _analytical(double px, double py, double pz, double t)
	{
		return sin(px*2.*freq*M_PI) * sin(py*2.*freq*M_PI) * sin(pz*2.*freq*M_PI) * exp(-4.*3*freq*freq*coeff*M_PI*M_PI*t);
	}
	
	double _analyticalRHS(double px, double py, double pz, double t)
	{
		return -freq*freq*4.*3.*M_PI*M_PI*_analytical(px,py,pz,t);
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
					{
						b(ix,iy,iz).tmpU = 0;
						b(ix,iy,iz).tmpV = 0;
						b(ix,iy,iz).tmpW = 0;
#ifdef _DENSITYDIFF_
						b(ix,iy,iz).tmp = 0;
#endif
					}
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
					{
						b(ix,iy,iz).u = b(ix,iy,iz).tmpU;
						b(ix,iy,iz).v = b(ix,iy,iz).tmpV;
						b(ix,iy,iz).w = b(ix,iy,iz).tmpW;
#ifdef _DENSITYDIFF_
						b(ix,iy,iz).rho = b(ix,iy,iz).tmp;
#endif
					}
		}
	}
	
	inline void diffuse(const double dt, const int stage)
	{
		const int N = vInfo.size();
		
#pragma omp parallel for schedule(static)
		for(int i=0; i<N; i++)
		{
			BlockInfo info = vInfo[i];
			FluidBlock& b = *(FluidBlock*)info.ptrBlock;
			const double prefactor = coeff * dt * ((stage==0)?.5:1);
			
			for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
				for(int iy=0; iy<FluidBlock::sizeY; ++iy)
					for(int ix=0; ix<FluidBlock::sizeX; ++ix)
					{
						double p[3];
						info.pos(p, ix, iy, iz);
						
						if (stage==0)
							b(ix,iy,iz).tmpU = b(ix,iy,iz).u + prefactor * _analyticalRHS(p[0],p[1],p[2],time);
						else if (stage==1)
							b(ix,iy,iz).tmpU = b(ix,iy,iz).u + prefactor * _analyticalRHS(p[0],p[1],p[2],time+dt*.5);
						b(ix,iy,iz).tmpV = b(ix,iy,iz).v;
						b(ix,iy,iz).tmpW = b(ix,iy,iz).w;
					}
		}
	}
	
public:
	CoordinatorDiffusionTimeTest(const double coeff, const double freq, FluidGrid * grid) : GenericCoordinator(grid), coeff(coeff), freq(freq), time(0)
	{
	}
	
	void operator()(const double dt)
	{
		check("diffusion - start");
		
		reset();
		diffuse(dt,0);
#ifdef _RK2_
		diffuse(dt,1);
#endif
		update();
		time+=dt;
		
		check("diffusion - end");
	}
	
	string getName()
	{
		return "DiffusionTimeTest";
	}
};
#endif
