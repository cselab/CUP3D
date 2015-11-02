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
    Real *uBody, *vBody;
	
	inline void reset()
	{
		const int N = vInfo.size();
		
#pragma omp parallel for schedule(static)
		for(int i=0; i<N; i++)
		{
			BlockInfo info = vInfo[i];
			FluidBlock& b = *(FluidBlock*)info.ptrBlock;
			
			for(int iy=0; iy<FluidBlock::sizeY; ++iy)
				for(int ix=0; ix<FluidBlock::sizeX; ++ix)
				{
					b(ix,iy).tmpU = 0;
					b(ix,iy).tmpV = 0;
#ifdef _DENSITYDIFF_
					b(ix,iy).tmp = 0;
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
			
			for(int iy=0; iy<FluidBlock::sizeY; ++iy)
				for(int ix=0; ix<FluidBlock::sizeX; ++ix)
				{
					b(ix,iy).u = b(ix,iy).tmpU;
					b(ix,iy).v = b(ix,iy).tmpV;
#ifdef _DENSITYDIFF_
					b(ix,iy).rho = b(ix,iy).tmp;
#endif
				}
		}
	 }
	
	inline void diffuse(const double dt, const int stage)
	{
		BlockInfo * ary = &vInfo.front();
		const int N = vInfo.size();
		
#pragma omp parallel
		{
			OperatorDiffusion kernel(dt, coeff, stage);
			//OperatorDiffusionHighOrder kernel(dt, coeff, stage);
			
            Lab mylab;
#ifdef _MOVING_FRAME_
            mylab.pDirichlet.u = 0;
            mylab.pDirichlet.v = *vBody;
#endif
			mylab.prepare(*grid, kernel.stencil_start, kernel.stencil_end, false);
			
#pragma omp for schedule(static)
			for (int i=0; i<N; i++)
			{
				mylab.load(ary[i], 0);
				kernel(mylab, ary[i], *(FluidBlock*)ary[i].ptrBlock);
			}
		}
	}
	
public:
	CoordinatorDiffusion(const double coeff, Real * uBody, Real * vBody, FluidGrid * grid) : GenericCoordinator(grid), coeff(coeff), uBody(uBody), vBody(vBody)
	{
	}
	
    CoordinatorDiffusion(const double coeff, FluidGrid * grid) : GenericCoordinator(grid), coeff(coeff), uBody(NULL), vBody(NULL)
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
    
    double _analytical(double px, double py, double t)
    {
        return sin(px*2.*freq*M_PI) * sin(py*2.*freq*M_PI) * exp(-4.*2*freq*freq*coeff*M_PI*M_PI*t);
    }
    
    double _analyticalRHS(double px, double py, double t)
    {
        return -freq*freq*4.*2.*M_PI*M_PI*_analytical(px,py,t);
    }
    
    inline void reset()
    {
        const int N = vInfo.size();
        
#pragma omp parallel for schedule(static)
        for(int i=0; i<N; i++)
        {
            BlockInfo info = vInfo[i];
            FluidBlock& b = *(FluidBlock*)info.ptrBlock;
            
            for(int iy=0; iy<FluidBlock::sizeY; ++iy)
                for(int ix=0; ix<FluidBlock::sizeX; ++ix)
                {
                    b(ix,iy).tmpU = 0;
                    b(ix,iy).tmpV = 0;
#ifdef _DENSITYDIFF_
                    b(ix,iy).tmp = 0;
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
            
            for(int iy=0; iy<FluidBlock::sizeY; ++iy)
                for(int ix=0; ix<FluidBlock::sizeX; ++ix)
                {
                    b(ix,iy).u = b(ix,iy).tmpU;
                    b(ix,iy).v = b(ix,iy).tmpV;
#ifdef _DENSITYDIFF_
                    b(ix,iy).rho = b(ix,iy).tmp;
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
            
            for(int iy=0; iy<FluidBlock::sizeY; ++iy)
                for(int ix=0; ix<FluidBlock::sizeX; ++ix)
                {
                    double p[3];
                    info.pos(p, ix, iy);
                    
                    if (stage==0)
                        b(ix,iy).tmpU = b(ix,iy).u + prefactor * _analyticalRHS(p[0],p[1],time);
                    else if (stage==1)
                        b(ix,iy).tmpU = b(ix,iy).u + prefactor * _analyticalRHS(p[0],p[1],time+dt*.5);
                    b(ix,iy).tmpV = b(ix,iy).v;
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
