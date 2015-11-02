//
//  OperatorAdvection.h
//  CubismUP_3D
//
//	Operates on
//		tmpU, tmpV
//
//  Created by Christian Conti on 1/7/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef CubismUP_3D_OperatorAdvection_h
#define CubismUP_3D_OperatorAdvection_h

#include <cmath>
#include "InterpolationKernels.h"
#include "GenericOperator.h"

template <typename RemeshingKernel>
class OperatorAdvection : public GenericLabOperator
{
private:
	double dt;
	
	template <typename Lab>
	Euler(Lab & lab, const BlockInfo& info, const int stage)
	{
		const double dh = info.h_gridpoint;
		const double invdh = 1./dh;
		
		const int support_start = RemeshingKernel::support_start;
		const int support_end   = RemeshingKernel::support_end;
		
		const int bx = info.index[0]*FluidBlock::sizeX;
		const int by = info.index[1]*FluidBlock::sizeY;
		
		for (int iy=support_start-1; iy<FluidBlock::sizeY+support_end; ++iy)
			for (int ix=support_start-1; ix<FluidBlock::sizeX+support_end; ++ix)
			{
				double p[2];
				info.pos(p,ix,iy);
			
				FluidElement& particle = lab(ix,iy);
			
				if (stage==0)
				{
#ifndef _RK2_
					particle.x = p[0] + dt * particle.u;
					particle.y = p[1] + dt * particle.v;
#else
					particle.x = p[0] + dt*.5 * particle.u;
					particle.y = p[1] + dt*.5 * particle.v;
#endif
				}
				else
				{
					particle.x = p[0] + dt * particle.tmpU;
					particle.y = p[1] + dt * particle.tmpV;
				}
			}
	}
	
	template <typename Lab>
	M2P(Lab & lab, const BlockInfo& info)
	{
		const double dh = info.h_gridpoint;
		const double invdh = 1./dh;
		
		const int support_start = RemeshingKernel::support_start;
		const int support_end   = RemeshingKernel::support_end;
		
		const int bx = info.index[0]*FluidBlock::sizeX;
		const int by = info.index[1]*FluidBlock::sizeY;
		
		for (int iy=support_start-1; iy<FluidBlock::sizeY+support_end; ++iy)
			for (int ix=support_start-1; ix<FluidBlock::sizeX+support_end; ++ix)
			{
				FluidElement& particle = lab(ix,iy);
				particle.tmpU = 0;
				particle.tmpV = 0;
			
				// nearest point with lower index
#ifndef _VERTEXCENTERED_
				double px = particle.x*invdh-.5;
				double py = particle.y*invdh-.5;
#else
				double px = particle.x*invdh;
				double py = particle.y*invdh;
#endif
				int fpx = (int)floor(px);
				int fpy = (int)floor(py);
			
				// compute weights
				double wx[RemeshingKernel::support], wy[RemeshingKernel::support];
				for (int i=support_start; i<support_end; i++)
				{
					wx[i-support_start] = RemeshingKernel::weight(px-(fpx+i));
					wy[i-support_start] = RemeshingKernel::weight(py-(fpy+i));
				}
			
				for (int j=support_start; j<support_end; j++)
					for (int i=support_start; i<support_end; i++)
					{
						const int lfpx = fpx+i - bx;
						const int lfpy = fpy+j - by;
						
						const double weight = wx[i-support_start] * wy[j-support_start];
				
						particle.tmpU += lab(lfpx,lfpy).u * weight;
						particle.tmpV += lab(lfpx,lfpy).v * weight;
					}
			}
	}
	
	template <typename Lab, typename BlockType>
	P2M(Lab & lab, const BlockInfo& info, BlockType& o)
	{
		const double dh = info.h_gridpoint;
		const double invdh = 1./dh;
		
		const int support_start = RemeshingKernel::support_start;
		const int support_end   = RemeshingKernel::support_end;
		
		const int bx = info.index[0]*FluidBlock::sizeX;
		const int by = info.index[1]*FluidBlock::sizeY;
		
		for (int iy=support_start-1; iy<FluidBlock::sizeY+support_end; ++iy)
		for (int ix=support_start-1; ix<FluidBlock::sizeX+support_end; ++ix)
		{
			double p[2];
			info.pos(p,ix,iy);
			
			FluidElement particle = lab(ix,iy);
			
			// P2M
			// nearest point with lower index
#ifndef _VERTEXCENTERED_
			const double px = particle.x*invdh-.5;
			const double py = particle.y*invdh-.5;
#else
			const double px = particle.x*invdh;
			const double py = particle.y*invdh;
#endif
			const int fpx = (int)floor(px);
			const int fpy = (int)floor(py);
			
			{
				// compute weights
				double wx[RemeshingKernel::support], wy[RemeshingKernel::support];
				for (int i=support_start; i<support_end; i++)
				{
					wx[i-support_start] = RemeshingKernel::weight(px-(fpx+i));
					wy[i-support_start] = RemeshingKernel::weight(py-(fpy+i));
				}
				
				// scatter only to elements within the block, elements outside the block are taken care by other blocks
				for (int j=support_start; j<support_end; j++)
				for (int i=support_start; i<support_end; i++)
				{
					if (fpx+i>=bx && fpx+i<bx+FluidBlock::sizeX &&
						fpy+j>=by && fpy+j<by+FluidBlock::sizeY)
					{
						const int lfpx = fpx+i - bx;
						const int lfpy = fpy+j - by;
						const double weight = wx[i-support_start] * wy[j-support_start];
						o(lfpx,lfpy).tmpU += weight * particle.u;
						o(lfpx,lfpy).tmpV += weight * particle.v;
#ifdef _MULTIPHASE_
						o(lfpx,lfpy).tmp += weight * particle.rho;
						//o(lfpx,lfpy).tmp += weight * particle.chi;
#endif
					}
				}
			}
			/*
			{
				// compute weights
				double wx[Hat::support], wy[Hat::support];
				for (int i=Hat::support_start; i<Hat::support_end; i++)
				{
					wx[i-Hat::support_start] = Hat::weight(px-(fpx+i));
					wy[i-Hat::support_start] = Hat::weight(py-(fpy+i));
				}
				
				// scatter only to elements within the block, elements outside the block are taken care by other blocks
				for (int j=Hat::support_start; j<Hat::support_end; j++)
				for (int i=Hat::support_start; i<Hat::support_end; i++)
				{
					if (fpx+i>=bx && fpx+i<bx+FluidBlock::sizeX &&
						fpy+j>=by && fpy+j<by+FluidBlock::sizeY)
					{
						const int lfpx = fpx+i - bx;
						const int lfpy = fpy+j - by;
						const double weight = wx[i-Hat::support_start] * wy[j-Hat::support_start];
#ifdef _MULTIPHASE_
						o(lfpx,lfpy).tmp += weight * lab(ix,iy).rho;
						//o(lfpx,lfpy).tmp += weight * lab(ix,iy).chi;
#endif
					}
				}
			}
			 */
		}
	}
	
public:
	OperatorAdvection(double dt) : dt(dt)
	{
		stencil_start[0] = -RemeshingKernel::support-2; // 1 for "CFL" and 1 for m2p
		stencil_start[1] = -RemeshingKernel::support-2;
		stencil_start[2] = 0;
		
		stencil_end[0] = RemeshingKernel::support+2;
		stencil_end[1] = RemeshingKernel::support+2;
		stencil_end[2] = 1;
	}
	~OperatorAdvection() {}
	
	template <typename Lab, typename BlockType>
	void operator()(Lab & lab, const BlockInfo& info, BlockType& o) const
	{
		//*
		Euler(lab, info, 0);
#ifndef _RK2_
		M2P(lab, info);
		Euler(lab, info, 1);
#endif
		P2M(lab, info, o);
		/*/
		const double dh = info.h_gridpoint;
		const double invdh = 1./dh;
		
		const int support_start = RemeshingKernel::support_start;
		const int support_end   = RemeshingKernel::support_end;
		
		const int bx = info.index[0]*FluidBlock::sizeX;
		const int by = info.index[1]*FluidBlock::sizeY;
		
		for (int iy=support_start-1; iy<FluidBlock::sizeY+support_end; ++iy)
			for (int ix=support_start-1; ix<FluidBlock::sizeX+support_end; ++ix)
			{
				double p[2];
				info.pos(p,ix,iy);
			
				FluidElement particle = lab(ix,iy);
			
				// RK2 midpoint - stage 1
				particle.x = p[0] + dt*.5 * particle.u;
				particle.y = p[1] + dt*.5 * particle.v;
			
				// M2P
				// nearest point with lower index
#ifndef _VERTEXCENTERED_
				double px = particle.x*invdh-.5;
				double py = particle.y*invdh-.5;
#else
				double px = particle.x*invdh;
				double py = particle.y*invdh;
#endif
				int fpx = (int)floor(px);
				int fpy = (int)floor(py);
			
				// compute weights
				double wx[RemeshingKernel::support], wy[RemeshingKernel::support];
				for (int i=support_start; i<support_end; i++)
				{
					wx[i-support_start] = RemeshingKernel::weight(px-(fpx+i));
					wy[i-support_start] = RemeshingKernel::weight(py-(fpy+i));
				}
			
				particle.u = 0;
				particle.v = 0;
				for (int j=support_start; j<support_end; j++)
					for (int i=support_start; i<support_end; i++)
					{
						const int lfpx = fpx+i - bx;
						const int lfpy = fpy+j - by;
						//if (lfpx < stencil_start[0] || lfpx >= bx+FluidBlock::sizeX+stencil_end[0]-1)
						//	cout << lfpx << " " << stencil_start[0] << " " << bx+FluidBlock::sizeX+stencil_end[0]-1;
						assert(lfpx >= stencil_start[0] && lfpx < bx+FluidBlock::sizeX+stencil_end[0]-1);
						//if (lfpy < stencil_start[1] || lfpy >= by+FluidBlock::sizeY+stencil_end[1]-1)
						//	cout << (p[1]-particle.y)*invdh << " " << py << " " << fpy << " " << j << " " << by << " " << lfpy << " " << stencil_start[1] << " " << by+FluidBlock::sizeY+stencil_end[1]-1;
						assert(lfpy >= stencil_start[1] && lfpy < by+FluidBlock::sizeY+stencil_end[1]-1);
						const double weight = wx[i-support_start] * wy[j-support_start];
				
						particle.u += lab(lfpx,lfpy).u * weight;
						particle.v += lab(lfpx,lfpy).v * weight;
					}
			
				// RK2 midpoint - stage 2
				particle.x = p[0] + dt * particle.u;
				particle.y = p[1] + dt * particle.v;
			
				// P2M
				// nearest point with lower index
#ifndef _VERTEXCENTERED_
				px = particle.x*invdh-.5;
				py = particle.y*invdh-.5;
#else
				px = particle.x*invdh;
				py = particle.y*invdh;
#endif
				fpx = (int)floor(px);
				fpy = (int)floor(py);
			
				// compute weights
				for (int i=support_start; i<support_end; i++)
				{
					wx[i-support_start] = RemeshingKernel::weight(px-(fpx+i));
					wy[i-support_start] = RemeshingKernel::weight(py-(fpy+i));
				}
			
				// scatter only to elements within the block, elements outside the block are taken care by other blocks
				for (int j=support_start; j<support_end; j++)
					for (int i=support_start; i<support_end; i++)
					{
						if (fpx+i>=bx && fpx+i<bx+FluidBlock::sizeX &&
						    fpy+j>=by && fpy+j<by+FluidBlock::sizeY)
						{
							const int lfpx = fpx+i - bx;
							const int lfpy = fpy+j - by;
							const double weight = wx[i-support_start] * wy[j-support_start];
					
							o(lfpx,lfpy).tmpU += weight * lab(ix,iy).u;
							o(lfpx,lfpy).tmpV += weight * lab(ix,iy).v;
#ifdef _MULTIPHASE_
							o(lfpx,lfpy).tmp += weight * lab(ix,iy).rho;
#endif // _MULTIPHASE_
						}
					}
		}
		 */
	}
};

template <typename RemeshingKernel>
class OperatorTransport : public GenericLabOperator
{
private:
	double dt;
	
    template <typename Lab>
    Euler(Lab & lab, const BlockInfo& info, const int stage)
    {
        const double dh = info.h_gridpoint;
        const double invdh = 1./dh;
        
        const int support_start = RemeshingKernel::support_start;
        const int support_end   = RemeshingKernel::support_end;
        
        const int bx = info.index[0]*FluidBlock::sizeX;
        const int by = info.index[1]*FluidBlock::sizeY;
        
        for (int iy=support_start-1; iy<FluidBlock::sizeY+support_end; ++iy)
            for (int ix=support_start-1; ix<FluidBlock::sizeX+support_end; ++ix)
            {
                double p[2];
                info.pos(p,ix,iy);
                
                FluidElement& particle = lab(ix,iy);
                
                if (stage==0)
                {
#ifndef _RK2_
                    particle.x = p[0] + dt * particle.u;
                    particle.y = p[1] + dt * particle.v;
#else
                    particle.x = p[0] + dt*.5 * particle.u;
                    particle.y = p[1] + dt*.5 * particle.v;
#endif
                }
                else
                {
                    particle.x = p[0] + dt * particle.tmpU;
                    particle.y = p[1] + dt * particle.tmpV;
                }
            }
    }
	
    template <typename Lab>
    M2P(Lab & lab, const BlockInfo& info)
    {
        const double dh = info.h_gridpoint;
        const double invdh = 1./dh;
        
        const int support_start = RemeshingKernel::support_start;
        const int support_end   = RemeshingKernel::support_end;
        
        const int bx = info.index[0]*FluidBlock::sizeX;
        const int by = info.index[1]*FluidBlock::sizeY;
        
        for(int iy=support_start-1; iy<FluidBlock::sizeY+support_end; ++iy)
            for(int ix=support_start-1; ix<FluidBlock::sizeX+support_end; ++ix)
            {
                FluidElement& particle = lab(ix,iy);
                particle.tmpU = 0;
                particle.tmpV = 0;
                
                // nearest point with lower index
#ifndef _VERTEXCENTERED_
                double px = particle.x*invdh-.5;
                double py = particle.y*invdh-.5;
#else
                double px = particle.x*invdh;
                double py = particle.y*invdh;
#endif
                int fpx = (int)floor(px);
                int fpy = (int)floor(py);
                
                // compute weights
                double wx[RemeshingKernel::support], wy[RemeshingKernel::support];
                for (int i=support_start; i<support_end; i++)
                {
                    wx[i-support_start] = RemeshingKernel::weight(px-(fpx+i));
                    wy[i-support_start] = RemeshingKernel::weight(py-(fpy+i));
                }
                
                for (int j=support_start; j<support_end; j++)
                    for (int i=support_start; i<support_end; i++)
                    {
                        const int lfpx = fpx+i - bx;
                        const int lfpy = fpy+j - by;
                        assert(lfpx >= stencil_start[0] && lfpx < FluidBlock::sizeX+stencil_end[0]-1);
                        assert(lfpy >= stencil_start[1] && lfpy < FluidBlock::sizeY+stencil_end[1]-1);
                        const double weight = wx[i-support_start] * wy[j-support_start];
                        
                        particle.tmpU += lab(lfpx,lfpy).u * weight;
                        particle.tmpV += lab(lfpx,lfpy).v * weight;
                    }
            }
    }
	
	template <typename Lab, typename BlockType>
    P2M(Lab & lab, const BlockInfo& info, BlockType& o)
    {
        const double dh = info.h_gridpoint;
        const double invdh = 1./dh;
        
        const int support_start = RemeshingKernel::support_start;
        const int support_end   = RemeshingKernel::support_end;
        
        const int bx = info.index[0]*FluidBlock::sizeX;
        const int by = info.index[1]*FluidBlock::sizeY;
        
        for (int iy=support_start-1; iy<FluidBlock::sizeY+support_end; ++iy)
            for (int ix=support_start-1; ix<FluidBlock::sizeX+support_end; ++ix)
            {
                double p[2];
                info.pos(p,ix,iy);
                
                FluidElement particle = lab(ix,iy);
                
                // P2M
                // nearest point with lower index
#ifndef _VERTEXCENTERED_
                const double px = particle.x*invdh-.5;
                const double py = particle.y*invdh-.5;
#else
                const double px = particle.x*invdh;
                const double py = particle.y*invdh;
#endif
                const int fpx = (int)floor(px);
                const int fpy = (int)floor(py);
                
                // compute weights
                double wx[RemeshingKernel::support], wy[RemeshingKernel::support];
                for (int i=support_start; i<support_end; i++)
                {
                    wx[i-support_start] = RemeshingKernel::weight(px-(fpx+i));
                    wy[i-support_start] = RemeshingKernel::weight(py-(fpy+i));
                }
                
                // scatter only to elements within the block, elements outside the block are taken care by other blocks
                for (int j=support_start; j<support_end; j++)
                    for (int i=support_start; i<support_end; i++)
                    {
                        if (fpx+i>=bx && fpx+i<bx+FluidBlock::sizeX &&
                            fpy+j>=by && fpy+j<by+FluidBlock::sizeY)
                        {
                            const int lfpx = fpx+i - bx;
                            const int lfpy = fpy+j - by;
                            const double weight = wx[i-support_start] * wy[j-support_start];
                            o(lfpx,lfpy).tmp += weight * lab(ix,iy).rho;
                        }
                    }
            }
	}
	
public:
	OperatorTransport(double dt) : dt(dt)
	{
		stencil_start[0] = -RemeshingKernel::support-2; // 1 for "CFL" and 1 for m2p
		stencil_start[1] = -RemeshingKernel::support-2;
		stencil_start[2] = 0;
		
		stencil_end[0] = RemeshingKernel::support+2;
		stencil_end[1] = RemeshingKernel::support+2;
		stencil_end[2] = 1;
	}
	
	~OperatorTransport() {}
	
	template <typename Lab, typename BlockType>
	void operator()(Lab & lab, const BlockInfo& info, BlockType& o) const
	{
		//*
		 Euler(lab, info, 0);
#ifdef _RK2_
		 M2P(lab, info);
		 Euler(lab, info, 1);
#endif
		 P2M(lab, info, o);
		 /*/
		const double dh = info.h_gridpoint;
		const double invdh = 1./dh;
		
		const int support_start = RemeshingKernel::support_start;
		const int support_end   = RemeshingKernel::support_end;
		
		const int bx = info.index[0]*FluidBlock::sizeX;
		const int by = info.index[1]*FluidBlock::sizeY;
		
		for(int iy=support_start-1; iy<FluidBlock::sizeY+support_end-2; ++iy)
		for(int ix=support_start-1; ix<FluidBlock::sizeX+support_end-2; ++ix)
		{
			double p[2];
			info.pos(p,ix,iy);
			
			FluidElement particle = lab(ix,iy);
			
			// RK2 midpoint - stage 1
			
			particle.x = p[0] + dt*.5 * particle.u;
			particle.y = p[1] + dt*.5 * particle.v;
			
			{
				// M2P
				// nearest point with lower index
#ifndef _VERTEXCENTERED_
				double px = particle.x*invdh-.5;
				double py = particle.y*invdh-.5;
#else
				double px = particle.x*invdh;
				double py = particle.y*invdh;
#endif
				int fpx = (int)floor(px);
				int fpy = (int)floor(py);
				
				// compute weights
				double wx[RemeshingKernel::support], wy[RemeshingKernel::support];
				for (int i=support_start; i<support_end; i++)
				{
					wx[i-support_start] = RemeshingKernel::weight(px-(fpx+i));
					wy[i-support_start] = RemeshingKernel::weight(py-(fpy+i));
				}
				
				for (int j=support_start; j<support_end; j++)
					for (int i=support_start; i<support_end; i++)
					{
						const int lfpx = fpx+i - bx;
						const int lfpy = fpy+j - by;
						assert(lfpx >= stencil_start[0] && lfpx < bx+FluidBlock::sizeX+stencil_end[0]-1);
						assert(lfpy >= stencil_start[1] && lfpy < by+FluidBlock::sizeY+stencil_end[1]-1);
						const double weight = wx[i-support_start] * wy[j-support_start];
				 
						particle.u += lab(lfpx,lfpy).u * weight;
						particle.v += lab(lfpx,lfpy).v * weight;
					}
			}
			
			// RK2 midpoint - stage 2
			particle.x = p[0] + dt * particle.u;
			particle.y = p[1] + dt * particle.v;
			
			{
				// P2M
				// nearest point with lower index
#ifndef _VERTEXCENTERED_
				double px = particle.x*invdh-.5;
				double py = particle.y*invdh-.5;
#else
				double px = particle.x*invdh;
				double py = particle.y*invdh;
#endif
				int fpx = (int)floor(px);
				int fpy = (int)floor(py);
				
				// compute weights
				double wx[RemeshingKernel::support], wy[RemeshingKernel::support];
				for (int i=support_start; i<support_end; i++)
				{
					wx[i-support_start] = RemeshingKernel::weight(px-(fpx+i));
					wy[i-support_start] = RemeshingKernel::weight(py-(fpy+i));
				}
				
				// scatter only to elements within the block, elements outside the block are taken care by other blocks
				for (int j=support_start; j<support_end; j++)
					for (int i=support_start; i<support_end; i++)
					{
						if (fpx+i>=bx && fpx+i<bx+FluidBlock::sizeX &&
							fpy+j>=by && fpy+j<by+FluidBlock::sizeY)
						{
							const int lfpx = fpx+i - bx;
							const int lfpy = fpy+j - by;
							const double weight = wx[i-support_start] * wy[j-support_start];
							o(lfpx,lfpy).tmp += weight * lab(ix,iy).rho;
						}
					}
			}
		}
		//*/
	}
};

template <typename RemeshingKernel>
class OperatorTransportTimeTest : public GenericLabOperator
{
private:
    double dt;
    double time;
    
    double _analyticalRHS(double px, double py, double t) const
    {
        return 8 * M_PI * cos((px-t) * 8. * M_PI);
    }
    
    template <typename Lab>
    Euler(Lab & lab, const BlockInfo& info, const int stage)
    {
        const double dh = info.h_gridpoint;
        const double invdh = 1./dh;
        
        const int support_start = RemeshingKernel::support_start;
        const int support_end   = RemeshingKernel::support_end;
        
        const int bx = info.index[0]*FluidBlock::sizeX;
        const int by = info.index[1]*FluidBlock::sizeY;
        
        for(int iy=support_start-1; iy<FluidBlock::sizeY+support_end; ++iy)
            for(int ix=support_start-1; ix<FluidBlock::sizeX+support_end; ++ix)
            {
                double p[2];
                info.pos(p,ix,iy);
                
                FluidElement& particle = lab(ix,iy);
                
                if (stage==0)
                {
#ifndef _RK2_
                    particle.x = p[0] + dt * 1;
                    particle.y = p[1] + dt * 1;
#else
                    particle.x = p[0] + dt*.5 * 1;
                    particle.y = p[1] + dt*.5 * 1;
#endif
                }
                else
                {
                    particle.x = p[0] + dt * 1;
                    particle.y = p[1] + dt * 1;
                }
            }
    }
    
    template <typename Lab>
    M2P(Lab & lab, const BlockInfo& info)
    {
        const double dh = info.h_gridpoint;
        const double invdh = 1./dh;
        
        const int support_start = RemeshingKernel::support_start;
        const int support_end   = RemeshingKernel::support_end;
        
        const int bx = info.index[0]*FluidBlock::sizeX;
        const int by = info.index[1]*FluidBlock::sizeY;
        
        for(int iy=support_start-1; iy<FluidBlock::sizeY+support_end; ++iy)
            for(int ix=support_start-1; ix<FluidBlock::sizeX+support_end; ++ix)
            {
                FluidElement& particle = lab(ix,iy);
                
                // nearest point with lower index
#ifndef _VERTEXCENTERED_
                double px = 2*(particle.x*invdh-.5)-1.;
                double py = 2*(particle.y*invdh-.5)-1.;
#else
                double px = 2*particle.x*invdh-1.;
                double py = 2*particle.y*invdh-1.;
#endif
                
                const Real r = sqrt(px*px + py*py);
                const Real invR = 1./r;
                
                particle.tmpU =  sin(py)*cos(r*M_PI/2)*invR;
                particle.tmpV = -sin(px)*cos(r*M_PI/2)*invR;
            }
    }
    
    template <typename Lab, typename BlockType>
    P2M(Lab & lab, const BlockInfo& info, BlockType& o)
    {
        const double dh = info.h_gridpoint;
        const double invdh = 1./dh;
        
        const int support_start = RemeshingKernel::support_start;
        const int support_end   = RemeshingKernel::support_end;
        
        const int bx = info.index[0]*FluidBlock::sizeX;
        const int by = info.index[1]*FluidBlock::sizeY;
        
        for(int iy=support_start-1; iy<FluidBlock::sizeY+support_end; ++iy)
            for(int ix=support_start-1; ix<FluidBlock::sizeX+support_end; ++ix)
            {
                double p[2];
                info.pos(p,ix,iy);
                
                FluidElement particle = lab(ix,iy);
                
                // P2M
                // nearest point with lower index
#ifndef _VERTEXCENTERED_
                const double px = particle.x*invdh-.5;
                const double py = particle.y*invdh-.5;
#else
                const double px = particle.x*invdh;
                const double py = particle.y*invdh;
#endif
                const int fpx = (int)floor(px);
                const int fpy = (int)floor(py);
                
                // compute weights
                double wx[RemeshingKernel::support], wy[RemeshingKernel::support];
                for (int i=support_start; i<support_end; i++)
                {
                    wx[i-support_start] = RemeshingKernel::weight(px-(fpx+i));
                    wy[i-support_start] = RemeshingKernel::weight(py-(fpy+i));
                }
                
                // scatter only to elements within the block, elements outside the block are taken care by other blocks
                for (int j=support_start; j<support_end; j++)
                    for (int i=support_start; i<support_end; i++)
                    {
                        if (fpx+i>=bx && fpx+i<bx+FluidBlock::sizeX &&
                            fpy+j>=by && fpy+j<by+FluidBlock::sizeY)
                        {
                            const int lfpx = fpx+i - bx;
                            const int lfpy = fpy+j - by;
                            const double weight = wx[i-support_start] * wy[j-support_start];
                            o(lfpx,lfpy).tmp += weight * lab(ix,iy).rho;
                        }
                    }
            }
    }
    
public:
    OperatorTransportTimeTest(double dt, double time) : dt(dt), time(time)
    {
        stencil_start[0] = -RemeshingKernel::support-2; // 1 for "CFL" and 1 for m2p
        stencil_start[1] = -RemeshingKernel::support-2;
        stencil_start[2] = 0;
        
        stencil_end[0] = RemeshingKernel::support+2;
        stencil_end[1] = RemeshingKernel::support+2;
        stencil_end[2] = 1;
    }
    
    ~OperatorTransportTimeTest() {}
    
    template <typename Lab, typename BlockType>
    void operator()(Lab & lab, const BlockInfo& info, BlockType& o) const
    {
        //*
        Euler(lab, info, 0);
#ifdef _RK2_
        M2P(lab, info);
        Euler(lab, info, 1);
#endif
        P2M(lab, info, o);
    }
};

class OperatorAdvectionFD : public GenericLabOperator
{
private:
	double dt;
	
public:
	OperatorAdvectionFD(double dt, const int stage) : dt(dt)
	{
		stencil_start[0] = -1;
		stencil_start[1] = -1;
		stencil_start[2] = 0;
		
		stencil_end[0] = 2;
		stencil_end[1] = 2;
		stencil_end[2] = 1;
	}
	~OperatorAdvectionFD() {}
	
	template <typename Lab, typename BlockType>
	void operator()(Lab & lab, const BlockInfo& info, BlockType& o) const
	{
		const double dh = info.h_gridpoint;
		const double invdh = dt*.5/dh;
		
		for (int iy=0; iy<FluidBlock::sizeY; ++iy)
			for (int ix=0; ix<FluidBlock::sizeX; ++ix)
			{
				o(ix,iy).tmpU = lab(ix,iy).u + lab(ix,iy).u * invdh * (lab(ix+1,iy).u - lab(ix-1,iy).u) + lab(ix,iy).v * invdh * (lab(ix,iy+1).u - lab(ix,iy-1).u);
				o(ix,iy).tmpV = lab(ix,iy).v + lab(ix,iy).u * invdh * (lab(ix+1,iy).v - lab(ix-1,iy).v) + lab(ix,iy).v * invdh * (lab(ix,iy+1).v - lab(ix,iy-1).v);
				o(ix,iy).tmp  = lab(ix,iy).rho + lab(ix,iy).u * invdh * (lab(ix+1,iy).rho - lab(ix-1,iy).rho) + lab(ix,iy).v * invdh * (lab(ix,iy+1).rho - lab(ix,iy-1).rho);
			}
	}
};

class OperatorAdvectionUpwind3rdOrder : public GenericLabOperator
{
private:
	double dt;
    const int stage;
    Real *uBody, *vBody;
	
public:
    OperatorAdvectionUpwind3rdOrder(double dt, Real * uBody, Real * vBody, const int stage) : dt(dt), uBody(uBody), vBody(vBody), stage(stage)
    {
        stencil_start[0] = -2;
        stencil_start[1] = -2;
        stencil_start[2] = 0;
        
        stencil_end[0] = 3;
        stencil_end[1] = 3;
        stencil_end[2] = 1;
    }
    
    OperatorAdvectionUpwind3rdOrder(double dt, const int stage) : dt(dt), uBody(NULL), vBody(NULL), stage(stage)
    {
        stencil_start[0] = -2;
        stencil_start[1] = -2;
        stencil_start[2] = 0;
        
        stencil_end[0] = 3;
        stencil_end[1] = 3;
        stencil_end[2] = 1;
    }
    
	~OperatorAdvectionUpwind3rdOrder() {}
	
	template <typename Lab, typename BlockType>
	void operator()(Lab & lab, const BlockInfo& info, BlockType& o) const
	{
#ifndef _RK2_
		const Real factor = -dt/(6.*info.h_gridpoint);
#else
		const Real factor = -dt*(stage==0 ? .5 : 1)/(6.*info.h_gridpoint);
#endif
		
		if (stage==0)
			for (int iy=0; iy<FluidBlock::sizeY; ++iy)
				for (int ix=0; ix<FluidBlock::sizeX; ++ix)
				{
					const Real dudx[2] = {  2*lab(ix+1,iy  ).u + 3*lab(ix  ,iy  ).u - 6*lab(ix-1,iy  ).u +   lab(ix-2,iy  ).u,
						                   -  lab(ix+2,iy  ).u + 6*lab(ix+1,iy  ).u - 3*lab(ix  ,iy  ).u - 2*lab(ix-1,iy  ).u};
			
					const Real dudy[2] = {  2*lab(ix  ,iy+1).u + 3*lab(ix  ,iy  ).u - 6*lab(ix  ,iy-1).u +   lab(ix  ,iy-2).u,
						                   -  lab(ix  ,iy+2).u + 6*lab(ix  ,iy+1).u - 3*lab(ix  ,iy  ).u - 2*lab(ix  ,iy-1).u};
			
					const Real dvdx[2] = {  2*lab(ix+1,iy  ).v + 3*lab(ix  ,iy  ).v - 6*lab(ix-1,iy  ).v +   lab(ix-2,iy  ).v,
						                   -  lab(ix+2,iy  ).v + 6*lab(ix+1,iy  ).v - 3*lab(ix  ,iy  ).v - 2*lab(ix-1,iy  ).v};
			
					const Real dvdy[2] = {  2*lab(ix  ,iy+1).v + 3*lab(ix  ,iy  ).v - 6*lab(ix  ,iy-1).v +   lab(ix  ,iy-2).v,
						                   -  lab(ix  ,iy+2).v + 6*lab(ix  ,iy+1).v - 3*lab(ix  ,iy  ).v - 2*lab(ix  ,iy-1).v};
			
					const Real u = o(ix,iy).u;
#ifndef _MOVING_FRAME_
                    const Real v = o(ix,iy).v;
#else
					const Real v = o(ix,iy).v - *vBody;
#endif
			
					o(ix,iy).tmpU = u + factor*(max(u,(Real)0) * dudx[0] + min(u,(Real)0) * dudx[1] +
												max(v,(Real)0) * dudy[0] + min(v,(Real)0) * dudy[1]);
					o(ix,iy).tmpV = o(ix,iy).v + factor*(max(u,(Real)0) * dvdx[0] + min(u,(Real)0) * dvdx[1] +
												max(v,(Real)0) * dvdy[0] + min(v,(Real)0) * dvdy[1]);
#ifndef _RK2_
#ifdef _MULTIPHASE_
					const Real drdx[2] = {  2*lab(ix+1,iy  ).rho + 3*lab(ix  ,iy  ).rho - 6*lab(ix-1,iy  ).rho +   lab(ix-2,iy  ).rho,
										   -  lab(ix+2,iy  ).rho + 6*lab(ix+1,iy  ).rho - 3*lab(ix  ,iy  ).rho - 2*lab(ix-1,iy  ).rho};
					
					const Real drdy[2] = {  2*lab(ix  ,iy+1).rho + 3*lab(ix  ,iy  ).rho - 6*lab(ix  ,iy-1).rho +   lab(ix  ,iy-2).rho,
									 	   -  lab(ix  ,iy+2).rho + 6*lab(ix  ,iy+1).rho - 3*lab(ix  ,iy  ).rho - 2*lab(ix  ,iy-1).rho};
					
					const Real r = o(ix,iy).rho;
					
					o(ix,iy).tmp  =r + factor*(max(u,(Real)0) * drdx[0] + min(u,(Real)0) * drdx[1] +
                                               max(v,(Real)0) * drdy[0] + min(v,(Real)0) * drdy[1]);
#endif // _MULTIPHASE_
#endif // _RK2_
			}
#ifdef _RK2_
		else if (stage==1)
			for (int iy=0; iy<FluidBlock::sizeY; ++iy)
				for (int ix=0; ix<FluidBlock::sizeX; ++ix)
				{
					const Real dudx[2] = {  2*lab(ix+1,iy  ).tmpU + 3*lab(ix  ,iy  ).tmpU - 6*lab(ix-1,iy  ).tmpU +   lab(ix-2,iy  ).tmpU,
									       -  lab(ix+2,iy  ).tmpU + 6*lab(ix+1,iy  ).tmpU - 3*lab(ix  ,iy  ).tmpU - 2*lab(ix-1,iy  ).tmpU};
			
					const Real dudy[2] = {  2*lab(ix  ,iy+1).tmpU + 3*lab(ix  ,iy  ).tmpU - 6*lab(ix  ,iy-1).tmpU +   lab(ix  ,iy-2).tmpU,
				                           -  lab(ix  ,iy+2).tmpU + 6*lab(ix  ,iy+1).tmpU - 3*lab(ix  ,iy  ).tmpU - 2*lab(ix  ,iy-1).tmpU};
			
					const Real dvdx[2] = {  2*lab(ix+1,iy  ).tmpV + 3*lab(ix  ,iy  ).tmpV - 6*lab(ix-1,iy  ).tmpV +   lab(ix-2,iy  ).tmpV,
										   -  lab(ix+2,iy  ).tmpV + 6*lab(ix+1,iy  ).tmpV - 3*lab(ix  ,iy  ).tmpV - 2*lab(ix-1,iy  ).tmpV};
			
					const Real dvdy[2] = {  2*lab(ix  ,iy+1).tmpV + 3*lab(ix  ,iy  ).tmpV - 6*lab(ix  ,iy-1).tmpV +   lab(ix  ,iy-2).tmpV,
							               -  lab(ix  ,iy+2).tmpV + 6*lab(ix  ,iy+1).tmpV - 3*lab(ix  ,iy  ).tmpV - 2*lab(ix  ,iy-1).tmpV};
			
					const Real drdx[2] = {  2*lab(ix+1,iy  ).tmp + 3*lab(ix  ,iy  ).tmp - 6*lab(ix-1,iy  ).tmp +   lab(ix-2,iy  ).tmp,
										   -  lab(ix+2,iy  ).tmp + 6*lab(ix+1,iy  ).tmp - 3*lab(ix  ,iy  ).tmp - 2*lab(ix-1,iy  ).tmp};
			
					const Real drdy[2] = {  2*lab(ix  ,iy+1).tmp + 3*lab(ix  ,iy  ).tmp - 6*lab(ix  ,iy-1).tmp +   lab(ix  ,iy-2).tmp,
										   -  lab(ix  ,iy+2).tmp + 6*lab(ix  ,iy+1).tmp - 3*lab(ix  ,iy  ).tmp - 2*lab(ix  ,iy-1).tmp};
			
                    const Real u = lab(ix,iy).tmpU;
#ifndef _MOVING_FRAME_
                    const Real v = lab(ix,iy).tmpV;
#else
                    const Real v = lab(ix,iy).tmpV - *vBody;
#endif
			
					o(ix,iy).tmpU = o(ix,iy).u + factor*(max(u,(Real)0) * dudx[0] + min(u,(Real)0) * dudx[1] +
                                                         max(v,(Real)0) * dudy[0] + min(v,(Real)0) * dudy[1]);
					o(ix,iy).tmpV = o(ix,iy).v + factor*(max(u,(Real)0) * dvdx[0] + min(u,(Real)0) * dvdx[1] +
										                 max(v,(Real)0) * dvdy[0] + min(v,(Real)0) * dvdy[1]);
#ifdef _MULTIPHASE_
                    const Real r = lab(ix,iy).rho;
					o(ix,iy).tmp = r + factor*(max(u,(Real)0) * drdx[0] + min(u,(Real)0) * drdx[1] +
                                               max(v,(Real)0) * drdy[0] + min(v,(Real)0) * drdy[1]);
#endif // _MULTIPHASE
                }
#endif // _RK2_
	}
};


class OperatorTransportUpwind3rdOrder : public GenericLabOperator
{
private:
	double dt;
	const int stage;
	
public:
	OperatorTransportUpwind3rdOrder(double dt, const int stage) : dt(dt), stage(stage)
	{
		stencil_start[0] = -2;
		stencil_start[1] = -2;
		stencil_start[2] = 0;
		
		stencil_end[0] = 3;
		stencil_end[1] = 3;
		stencil_end[2] = 1;
	}
    ~OperatorTransportUpwind3rdOrder() {}
    
    template <typename Lab, typename BlockType>
    void operator()(Lab & lab, const BlockInfo& info, BlockType& o) const
    {
#ifndef _RK2_
        const Real factor = -dt/(6.*info.h_gridpoint);
#else
        const Real factor = -dt*(stage==0 ? .5 : 1)/(6.*info.h_gridpoint);
#endif
        
        if (stage==0)
            for (int iy=0; iy<FluidBlock::sizeY; ++iy)
                for (int ix=0; ix<FluidBlock::sizeX; ++ix)
                {
                    double p[3];
                    info.pos(p, ix, iy);
                    
                    const Real drdx[2] = {  2*lab(ix+1,iy  ).rho + 3*lab(ix  ,iy  ).rho - 6*lab(ix-1,iy  ).rho +   lab(ix-2,iy  ).rho,
                                           -  lab(ix+2,iy  ).rho + 6*lab(ix+1,iy  ).rho - 3*lab(ix  ,iy  ).rho - 2*lab(ix-1,iy  ).rho};
                    
                    const Real drdy[2] = {  2*lab(ix  ,iy+1).rho + 3*lab(ix  ,iy  ).rho - 6*lab(ix  ,iy-1).rho +   lab(ix  ,iy-2).rho,
                                           -  lab(ix  ,iy+2).rho + 6*lab(ix  ,iy+1).rho - 3*lab(ix  ,iy  ).rho - 2*lab(ix  ,iy-1).rho};
                    
                    const Real u = o(ix,iy).u;
                    const Real v = o(ix,iy).v;
                    const Real r = o(ix,iy).rho;
#ifndef _RK2_
                    o(ix,iy).tmp  = r + factor*(max(u,(Real)0) * drdx[0] + min(u,(Real)0) * drdx[1] +
                                                max(v,(Real)0) * drdy[0] + min(v,(Real)0) * drdy[1]);
#else
                    const Real dudx[2] = {  2*lab(ix+1,iy  ).u + 3*lab(ix  ,iy  ).u - 6*lab(ix-1,iy  ).u +   lab(ix-2,iy  ).u,
                                           -  lab(ix+2,iy  ).u + 6*lab(ix+1,iy  ).u - 3*lab(ix  ,iy  ).u - 2*lab(ix-1,iy  ).u};
                    
                    const Real dudy[2] = {  2*lab(ix  ,iy+1).u + 3*lab(ix  ,iy  ).u - 6*lab(ix  ,iy-1).u +   lab(ix  ,iy-2).u,
                                           -  lab(ix  ,iy+2).u + 6*lab(ix  ,iy+1).u - 3*lab(ix  ,iy  ).u - 2*lab(ix  ,iy-1).u};
                    
                    const Real dvdx[2] = {  2*lab(ix+1,iy  ).v + 3*lab(ix  ,iy  ).v - 6*lab(ix-1,iy  ).v +   lab(ix-2,iy  ).v,
                                           -  lab(ix+2,iy  ).v + 6*lab(ix+1,iy  ).v - 3*lab(ix  ,iy  ).v - 2*lab(ix-1,iy  ).v};
                    
                    const Real dvdy[2] = {  2*lab(ix  ,iy+1).v + 3*lab(ix  ,iy  ).v - 6*lab(ix  ,iy-1).v +   lab(ix  ,iy-2).v,
                                           -  lab(ix  ,iy+2).v + 6*lab(ix  ,iy+1).v - 3*lab(ix  ,iy  ).v - 2*lab(ix  ,iy-1).v};
                    
                    o(ix,iy).tmpU = u + factor*(max(u,(Real)0) * dudx[0] + min(u,(Real)0) * dudx[1] +
                                                max(v,(Real)0) * dudy[0] + min(v,(Real)0) * dudy[1]);
                    o(ix,iy).tmpV = v + factor*(max(u,(Real)0) * dvdx[0] + min(u,(Real)0) * dvdx[1] +
                                                max(v,(Real)0) * dvdy[0] + min(v,(Real)0) * dvdy[1]);
                    o(ix,iy).tmp  = r + factor*(max(u,(Real)0) * drdx[0] + min(u,(Real)0) * drdx[1] +
                                                max(v,(Real)0) * drdy[0] + min(v,(Real)0) * drdy[1]);
#endif
                }
#ifdef _RK2_
        else if (stage==1)
            for (int iy=0; iy<FluidBlock::sizeY; ++iy)
                for (int ix=0; ix<FluidBlock::sizeX; ++ix)
                {
                    double p[3];
                    info.pos(p, ix, iy);
                    
                    const Real drdx[2] = {  2*lab(ix+1,iy  ).tmp + 3*lab(ix  ,iy  ).tmp - 6*lab(ix-1,iy  ).tmp +   lab(ix-2,iy  ).tmp,
                                           -  lab(ix+2,iy  ).tmp + 6*lab(ix+1,iy  ).tmp - 3*lab(ix  ,iy  ).tmp - 2*lab(ix-1,iy  ).tmp};
                    
                    const Real drdy[2] = {  2*lab(ix  ,iy+1).tmp + 3*lab(ix  ,iy  ).tmp - 6*lab(ix  ,iy-1).tmp +   lab(ix  ,iy-2).tmp,
                                           -  lab(ix  ,iy+2).tmp + 6*lab(ix  ,iy+1).tmp - 3*lab(ix  ,iy  ).tmp - 2*lab(ix  ,iy-1).tmp};
                    
                    const Real u = lab(ix,iy).tmpU;
                    const Real v = lab(ix,iy).tmpV;
                    const Real r = lab(ix,iy).rho;
                    
                    o(ix,iy).tmp  = r + factor*(max(u,(Real)0) * drdx[0] + min(u,(Real)0) * drdx[1] +
                                                max(v,(Real)0) * drdy[0] + min(v,(Real)0) * drdy[1]);
                }
#endif
    }
};


class OperatorTransportTimeTestUpwind3rdOrder : public GenericLabOperator
{
private:
    double dt;
    double time;
    const int stage;
    
    double _analyticalRHS(double px, double py, double t) const
    {
        return 8 * M_PI * cos((px-t) * 8. * M_PI);
    }
    
public:
    OperatorTransportTimeTestUpwind3rdOrder(double dt, double time, const int stage) : dt(dt), time(time), stage(stage)
    {
        stencil_start[0] = -2;
        stencil_start[1] = -2;
        stencil_start[2] = 0;
        
        stencil_end[0] = 3;
        stencil_end[1] = 3;
        stencil_end[2] = 1;
    }
    ~OperatorTransportTimeTestUpwind3rdOrder() {}
    
    template <typename Lab, typename BlockType>
    void operator()(Lab & lab, const BlockInfo& info, BlockType& o) const
    {
#ifndef _RK2_
        const double factor = -dt;
#else
        const double factor = -dt * ((stage==0)?.5:1);
#endif
        
        if (stage==0)
            for (int iy=0; iy<FluidBlock::sizeY; ++iy)
                for (int ix=0; ix<FluidBlock::sizeX; ++ix)
                {
                    double p[3];
                    info.pos(p, ix, iy);
                    
                    const Real drdx[2] = { _analyticalRHS(p[0],p[1],time),_analyticalRHS(p[0],p[1],time) };
                    const Real drdy[2] = { 0,0 };
                    
                    const Real u = o(ix,iy).u;
                    const Real v = o(ix,iy).v;
                    const Real r = o(ix,iy).rho;
#ifndef _RK2_
                    o(ix,iy).tmp  = r + factor*(max(u,(Real)0) * drdx[0] + min(u,(Real)0) * drdx[1] +
                                                max(v,(Real)0) * drdy[0] + min(v,(Real)0) * drdy[1]);
#else
                    const Real dudx[2] = {  2*lab(ix+1,iy  ).u + 3*lab(ix  ,iy  ).u - 6*lab(ix-1,iy  ).u +   lab(ix-2,iy  ).u,
                                           -  lab(ix+2,iy  ).u + 6*lab(ix+1,iy  ).u - 3*lab(ix  ,iy  ).u - 2*lab(ix-1,iy  ).u};
                    
                    const Real dudy[2] = {  2*lab(ix  ,iy+1).u + 3*lab(ix  ,iy  ).u - 6*lab(ix  ,iy-1).u +   lab(ix  ,iy-2).u,
                                           -  lab(ix  ,iy+2).u + 6*lab(ix  ,iy+1).u - 3*lab(ix  ,iy  ).u - 2*lab(ix  ,iy-1).u};
                    
                    const Real dvdx[2] = {  2*lab(ix+1,iy  ).v + 3*lab(ix  ,iy  ).v - 6*lab(ix-1,iy  ).v +   lab(ix-2,iy  ).v,
                                           -  lab(ix+2,iy  ).v + 6*lab(ix+1,iy  ).v - 3*lab(ix  ,iy  ).v - 2*lab(ix-1,iy  ).v};
                    
                    const Real dvdy[2] = {  2*lab(ix  ,iy+1).v + 3*lab(ix  ,iy  ).v - 6*lab(ix  ,iy-1).v +   lab(ix  ,iy-2).v,
                                           -  lab(ix  ,iy+2).v + 6*lab(ix  ,iy+1).v - 3*lab(ix  ,iy  ).v - 2*lab(ix  ,iy-1).v};
                    
                    o(ix,iy).tmpU = u + factor*(max(u,(Real)0) * dudx[0] + min(u,(Real)0) * dudx[1] +
                                                max(v,(Real)0) * dudy[0] + min(v,(Real)0) * dudy[1]);
                    o(ix,iy).tmpV = v + factor*(max(u,(Real)0) * dvdx[0] + min(u,(Real)0) * dvdx[1] +
                                                max(v,(Real)0) * dvdy[0] + min(v,(Real)0) * dvdy[1]);
                    o(ix,iy).tmp  = r + factor*(max(u,(Real)0) * drdx[0] + min(u,(Real)0) * drdx[1] +
                                                max(v,(Real)0) * drdy[0] + min(v,(Real)0) * drdy[1]);
#endif
                }
#ifdef _RK2_
        else if (stage==1)
            for (int iy=0; iy<FluidBlock::sizeY; ++iy)
                for (int ix=0; ix<FluidBlock::sizeX; ++ix)
                {
                    double p[3];
                    info.pos(p, ix, iy);
                    
                    const Real drdx[2] = { _analyticalRHS(p[0],p[1],time+dt*.5),_analyticalRHS(p[0],p[1],time+dt*.5) };
                    const Real drdy[2] = { 0,0 };
                    
                    const Real u = lab(ix,iy).tmpU;
                    const Real v = lab(ix,iy).tmpV;
                    const Real r = lab(ix,iy).rho;
                    
                    o(ix,iy).tmp  = r + factor*(max(u,(Real)0) * drdx[0] + min(u,(Real)0) * drdx[1] +
                                                max(v,(Real)0) * drdy[0] + min(v,(Real)0) * drdy[1]);
                }
#endif
    }
};

#endif
