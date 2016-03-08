//
//  OperatorMovingPressure.h
//  CubismUP_3D
//
//  Created by Christian Conti on 3/8/16.
//  Copyright © 2016 ETHZ. All rights reserved.
//

#ifndef OperatorMovingPressure_h
#define OperatorMovingPressure_h

#include "GenericOperator.h"
#include "InterpolationKernels.h"

class OperatorMovingPressure : public GenericLabOperator
{
private:
	double dt;
	Real uBody, vBody, wBody;
	
public:
	OperatorMovingPressure(double dt, Real uBody, Real vBody, Real wBody) : dt(dt), uBody(uBody), vBody(vBody), wBody(wBody)
	{
		// Mp4 stencil+backward step with CFL<1
		stencil = StencilInfo(-2,-2,-2, 3,3,3, false, 1, 6);
		
		stencil_start[0] = -2;
		stencil_start[1] = -2;
		stencil_start[2] = -2;
		stencil_end[0] = 3;
		stencil_end[1] = 3;
		stencil_end[2] = 3;
	}
	
	~OperatorMovingPressure() {}
	
	template <typename Lab, typename BlockType>
	void operator()(Lab & lab, const BlockInfo& info, BlockType& o) const
	{
		const double dh = info.h_gridpoint;
		const double invdh = 1./dh;
		
		const int bx = info.index[0]*FluidBlock::sizeX;
		const int by = info.index[1]*FluidBlock::sizeY;
		const int bz = info.index[2]*FluidBlock::sizeZ;
		
		for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
			for(int iy=0; iy<FluidBlock::sizeY; ++iy)
				for(int ix=0; ix<FluidBlock::sizeX; ++ix)
				{
					double p[3];
					info.pos(p,ix,iy,iz);
					
					const FluidElement& particle = lab(ix,iy,iz);
					o(ix,iy,iz).divU = 0;
					
					p[0] -= dt*uBody;
					p[1] -= dt*vBody;
					p[2] -= dt*wBody;
					
					
					// P2M
					// nearest point with lower index
#ifndef _VERTEXCENTERED_
					const double px = p[0]*invdh-.5;
					const double py = p[1]*invdh-.5;
					const double pz = p[2]*invdh-.5;
#else
					const double px = p[0]*invdh;
					const double py = p[1]*invdh;
					const double pz = p[2]*invdh;
#endif
					const int fpx = (int)floor(px);
					const int fpy = (int)floor(py);
					const int fpz = (int)floor(pz);
					
					{
						// compute weights
						double wx[Mp4::support], wy[Mp4::support], wz[Mp4::support];
						for (int i=Mp4::support_start; i<Mp4::support_end; i++)
						{
							wx[i-Mp4::support_start] = Mp4::weight(px-(fpx+i));
							wy[i-Mp4::support_start] = Mp4::weight(py-(fpy+i));
							wz[i-Mp4::support_start] = Mp4::weight(pz-(fpz+i));
						}
						
						// scatter only to elements within the block, elements outside the block are taken care by other blocks
						for (int k=Mp4::support_start; k<Mp4::support_end; k++)
							for (int j=Mp4::support_start; j<Mp4::support_end; j++)
								for (int i=Mp4::support_start; i<Mp4::support_end; i++)
								{
									const int lfpx = fpx+i - bx;
									const int lfpy = fpy+j - by;
									const int lfpz = fpz+k - bz;
									
									const double weight = wx[i-Mp4::support_start] * wy[j-Mp4::support_start] * wz[k-Mp4::support_start];
									
									o(ix,iy,iz).divU += lab(lfpx,lfpy,lfpz).pOld * weight;
								}
					}
				}
	}
};

#endif /* OperatorMovingPressure_h */
