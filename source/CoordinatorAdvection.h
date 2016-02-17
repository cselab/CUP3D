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
	
#ifdef _ISPC_
	enum {
		_BSX2_ = _BSX_+4,
		_BSY2_ = _BSY_+4,
		_BSZ2_ = _BSZ_+4,
		SLICESIZE = _BSX_*_BSY_,
		SLICESIZE2 = (_BSX_+4)*(_BSY_+4),
		SIZE = _BSX_*_BSY_*_BSZ_,
		SIZE2 = (_BSX_+4)*(_BSY_+4)*(_BSZ_+4)
	};
	
	Real * uLab;
	Real * vLab;
	Real * wLab;
	Real * u;
	Real * v;
	Real * w;
	Real * uTmp;
	Real * vTmp;
	Real * wTmp;
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
#ifndef _ISPC_
		OperatorAdvectionUpwind3rdOrder kernel(dt,uBody,vBody,wBody,stage);
#else // _ISPC_
#ifdef _MULTIPHASE_
#warning ISPC with MULTIPHASE unsupported
		cout << "ISPC advection with MULTIPHASE unsupported yet!\n";
		abort();
#else // _MULTIPHASE_
#ifndef _RK2_
		OperatorAdvectionUpwind3rdOrderISPC kernel(dt,uBody,vBody,wBody,u,v,w,uTmp,vTmp,wTmp,uLab,vLab,wLab,stage);
#else // _RK2_
#warning ISPC with RK2 unsupported
		cout << "ISPC advection with RK2 unsupported yet!\n";
		abort();
#endif // _RK2_
#endif // _MULTIPHASE_
#endif // _ISPC_
		compute(kernel);
	}
	
public:
#ifndef _MULTIPHASE_
	CoordinatorAdvection(Real * uBody, Real * vBody, Real * wBody, FluidGridMPI * grid) : GenericCoordinator(grid), uBody(uBody), vBody(vBody), wBody(wBody)
#else
	CoordinatorAdvection(Real * uBody, Real * vBody, Real * wBody, FluidGridMPI * grid, Real rhoS) : GenericCoordinator(grid), uBody(uBody), vBody(vBody), wBody(wBody), rhoS(rhoS)
#endif
	{
#ifdef _ISPC_
		uLab = new Real[SIZE2*NTHREADS];
		vLab = new Real[SIZE2*NTHREADS];
		wLab = new Real[SIZE2*NTHREADS];
		u = new Real[SIZE*NTHREADS];
		v = new Real[SIZE*NTHREADS];
		w = new Real[SIZE*NTHREADS];
		uTmp = new Real[SIZE*NTHREADS];
		vTmp = new Real[SIZE*NTHREADS];
		wTmp = new Real[SIZE*NTHREADS];
#endif
	}
	
#ifndef _MULTIPHASE_
	CoordinatorAdvection(FluidGridMPI * grid) : GenericCoordinator(grid), uBody(NULL), vBody(NULL)
#else
	CoordinatorAdvection(FluidGridMPI * grid, Real rhoS) : GenericCoordinator(grid), uBody(NULL), vBody(NULL), wBody(NULL), rhoS(rhoS)
#endif
	{
#ifdef _ISPC_
		uLab = new Real[SIZE2*NTHREADS];
		vLab = new Real[SIZE2*NTHREADS];
		wLab = new Real[SIZE2*NTHREADS];
		u = new Real[SIZE*NTHREADS];
		v = new Real[SIZE*NTHREADS];
		w = new Real[SIZE*NTHREADS];
		uTmp = new Real[SIZE*NTHREADS];
		vTmp = new Real[SIZE*NTHREADS];
		wTmp = new Real[SIZE*NTHREADS];
#endif
	}
	
	~CoordinatorAdvection()
	{
#ifdef _ISPC_
		delete [] uLab;
		delete [] vLab;
		delete [] wLab;
		delete [] u;
		delete [] v;
		delete [] w;
		delete [] uTmp;
		delete [] vTmp;
		delete [] wTmp;
#endif
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
	CoordinatorTransport(FluidGridMPI * grid) : GenericCoordinator(grid)
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
	CoordinatorTransportTimeTest(FluidGridMPI * grid) : GenericCoordinator(grid)
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
