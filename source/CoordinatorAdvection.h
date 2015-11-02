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
#ifdef _MULTIPHASE_
					b(ix,iy,iz).tmp = 0;
#endif // _MULTIPHASE_
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
#ifdef _MULTIPHASE_
						//b(ix,iy).chi = b(ix,iy).tmp;
						//b(ix,iy).rho = b(ix,iy).chi * rhoS + (1-b(ix,iy).chi);
						
						// threshold density
#ifdef _PARTICLES_
						Real density = min(max(b(ix,iy,iz).tmp,min((Real)1.,rhoS)),max((Real)1.,rhoS));
						b(ix,iy,iz).rho = density;
#else // _PARTICLES_
						//b(ix,iy).rho = b(ix,iy).tmp;
						Real density = min(max(b(ix,iy,iz).tmp,min((Real)1.,rhoS)),max((Real)1.,rhoS));
						b(ix,iy,iz).rho = density;
#endif // _PARTICLES_
#endif // _MULTIPHASE_
					}
			}
	}
	
	inline void advect(const double dt)
	{
		BlockInfo * ary = &vInfo.front();
		const int N = vInfo.size();
		
#pragma omp parallel
		{
#ifndef _PARTICLES_
			OperatorAdvectionUpwind3rdOrder kernel(dt,uBody,vBody,wBody,0);
			//OperatorAdvectionFD kernel(dt);
			
            Lab mylab;
#ifdef _MOVING_FRAME_
			mylab.pDirichlet.u = 0;
			mylab.pDirichlet.v = *vBody;
			mylab.pDirichlet.w = 0;
#endif
			mylab.prepare(*grid, kernel.stencil_start, kernel.stencil_end, false);
#else // _PARTICLES_
			//OperatorAdvection<Hat> kernel(dt);
			//OperatorAdvection<Lambda2> kernel(dt);
			//OperatorAdvection<Mp4> kernel(dt);
			OperatorAdvection<Ms6> kernel(dt);
			
            Lab mylab;
#ifdef _MOVING_FRAME_
			mylab.pDirichlet.u = 0;
			mylab.pDirichlet.v = *vBody;
			mylab.pDirichlet.w = 0;
#endif
			mylab.prepare(*grid, kernel.stencil_start, kernel.stencil_end, true);
#endif // _PARTICLES_
			
#pragma omp for schedule(static)
			for (int i=0; i<N; i++)
			{
				mylab.load(ary[i], 0);
				
				kernel(mylab, ary[i], *(FluidBlock*)ary[i].ptrBlock);
			}
		}
        
#ifndef _PARTICLES_
#ifdef _RK2_
#pragma omp parallel
		{
			// this is wrong - using -u instead of u?
			OperatorAdvectionUpwind3rdOrder kernel(dt,uBody,vBody,wBody,1);
			//OperatorAdvectionFD kernel(dt);
			
            Lab mylab;
#ifdef _MOVING_FRAME_
            mylab.pDirichlet.u = 0;
			mylab.pDirichlet.v = *vBody;
			mylab.pDirichlet.w = 0;
#endif
			mylab.prepare(*grid, kernel.stencil_start, kernel.stencil_end, false);
			
#pragma omp for schedule(static)
			for (int i=0; i<N; i++)
			{
				mylab.load(ary[i], 0);
				
				kernel(mylab, ary[i], *(FluidBlock*)ary[i].ptrBlock);
			}
		}
#endif // _RK2_
#endif // _PARTICLES_
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

		reset();
		advect(dt);
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
#ifndef _PARTICLES_
            OperatorTransportUpwind3rdOrder kernel(dt,0);
            
            Lab mylab;
            mylab.prepare(*grid, kernel.stencil_start, kernel.stencil_end, false);
#else
            //OperatorTransport<Hat> kernel(dt);
            OperatorTransport<Mp4> kernel(dt);
            //OperatorTransport<Ms6> kernel(dt);
            
            Lab mylab;
            mylab.prepare(*grid, kernel.stencil_start, kernel.stencil_end, true);
#endif
            
#pragma omp for schedule(static)
            for (int i=0; i<N; i++)
            {
                mylab.load(ary[i], 0);
                
                kernel(mylab, ary[i], *(FluidBlock*)ary[i].ptrBlock);
            }
        }
        
#ifndef _PARTICLES_
#ifdef _RK2_
#pragma omp parallel
        {
            OperatorTransportUpwind3rdOrder kernel(dt,1);
            
            Lab mylab;
            mylab.prepare(*grid, kernel.stencil_start, kernel.stencil_end, true);
            
#pragma omp for schedule(static)
            for (int i=0; i<N; i++)
            {
                mylab.load(ary[i], 0);
                
                kernel(mylab, ary[i], *(FluidBlock*)ary[i].ptrBlock);
            }
        }
#endif // _RK2_
#endif // _PARTICLES_
        
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
        
        /*
#pragma omp parallel for schedule(static)
        for(int i=0; i<N; i++)
        {
            BlockInfo info = vInfo[i];
            FluidBlock& b = *(FluidBlock*)info.ptrBlock;
#ifndef _RK2_
            const double prefactor = dt;
#else
            const double prefactor = dt * ((stage==0)?.5:1);
#endif
            
            for(int iy=0; iy<FluidBlock::sizeY; ++iy)
                for(int ix=0; ix<FluidBlock::sizeX; ++ix)
                {
                    double p[3];
                    info.pos(p, ix, iy);
                    
                    if (stage==0)
                        b(ix,iy).tmp = b(ix,iy).rho + prefactor * _analyticalRHS(p[0], p[1], time);
                    else if (stage==1)
                        b(ix,iy).tmp = b(ix,iy).rho + prefactor * _analyticalRHS(p[0], p[1], time+dt*.5);
                }
        }
        
         /*/
        if (stage==0)
#pragma omp parallel
        {
#ifndef _PARTICLES_
            OperatorTransportUpwind3rdOrder kernel(dt,0);
            //OperatorTransportTimeTestUpwind3rdOrder kernel(dt,time,0);
            
            Lab mylab;
            mylab.prepare(*grid, kernel.stencil_start, kernel.stencil_end, false);
#else
            //OperatorTransportTimeTest<Mp4> kernel(dt,time);
            //OperatorTransportTimeTest<Ms6> kernel(dt,time);
            OperatorTransport<Ms6> kernel(dt);
            
            Lab mylab;
            mylab.prepare(*grid, kernel.stencil_start, kernel.stencil_end, true);
#endif
            
#pragma omp for schedule(static)
            for (int i=0; i<N; i++)
            {
                mylab.load(ary[i], 0);
                
                kernel(mylab, ary[i], *(FluidBlock*)ary[i].ptrBlock);
            }
        }
#ifndef _PARTICLES_
#ifdef _RK2_
        else if (stage==1)
#pragma omp parallel
        {
            OperatorTransportUpwind3rdOrder kernel(dt,1);
            //OperatorTransportTimeTestUpwind3rdOrder kernel(dt,time,1);
            
            Lab mylab;
            mylab.prepare(*grid, kernel.stencil_start, kernel.stencil_end, false);
            
#pragma omp for schedule(static)
            for (int i=0; i<N; i++)
            {
                mylab.load(ary[i], 0);
                
                kernel(mylab, ary[i], *(FluidBlock*)ary[i].ptrBlock);
            }
        }
#endif
#endif
        //*/
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
#ifndef _PARTICLES_
#ifdef _RK2_
        advect(dt,1);
#endif
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
