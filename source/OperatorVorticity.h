//
//  OperatorVorticity.h
//  CubismUP_3D
//
//	Operates on
//		tmp
//
//  Created by Christian Conti on 1/12/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef CubismUP_3D_OperatorVorticity_h
#define CubismUP_3D_OperatorVorticity_h

#include "Layer.h"
#include "GenericOperator.h"

class OperatorVorticity : public GenericLabOperator
{
private:
	Layer & vorticity;
	
public:
	OperatorVorticity(Layer & vorticity) : vorticity(vorticity)
	{
		stencil_start[0] = -1;
		stencil_start[1] = -1;
		stencil_start[2] = 0;
		stencil_end[0] = 2;
		stencil_end[1] = 2;
		stencil_end[2] = 1;
	}
	
	~OperatorVorticity() {}
	
	template <typename Lab, typename BlockType>
	void operator()(Lab & lab, const BlockInfo& info, BlockType& o) const
	{
		const Real factor = 0.5/info.h_gridpoint;
		const int bx = info.index[0]*FluidBlock::sizeX;
		const int by = info.index[1]*FluidBlock::sizeY;
		
		for(int iy=0; iy<FluidBlock::sizeY; ++iy)
			for(int ix=0; ix<FluidBlock::sizeX; ++ix)
			{
				const Real vW = lab(ix-1,iy).v;
				const Real vE = lab(ix+1,iy).v;
				const Real uS = lab(ix,iy-1).u;
				const Real uN = lab(ix,iy+1).u;
				
				o(ix,iy).tmp = factor * ((vE-vW) - (uN-uS));
				vorticity(bx + ix, by + iy) = factor * ((vE-vW) - (uN-uS));
			}
	}
};

#endif
