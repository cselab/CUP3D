//
//  OperatorVorticity.h
//  CubismUP_3D
//
//  Created by Christian Conti on 12/8/15.
//  Copyright © 2015 ETHZ. All rights reserved.
//

#ifndef OperatorVorticity_h
#define OperatorVorticity_h

#include "GenericOperator.h"

class OperatorVorticity : public GenericLabOperator
{
private:
	
public:
	OperatorVorticity()
	{
		stencil_start[0] = -1;
		stencil_start[1] = -1;
		stencil_start[2] = -1;
		stencil_end[0] = 2;
		stencil_end[1] = 2;
		stencil_end[2] = 2;
	}
	
	~OperatorVorticity() {}
	
	template <typename Lab, typename BlockType>
	void operator()(Lab & lab, const BlockInfo& info, BlockType& o) const
	{
		const double inv2h = .5 / info.h_gridpoint;
		
		for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
		for(int iy=0; iy<FluidBlock::sizeY; ++iy)
		for(int ix=0; ix<FluidBlock::sizeX; ++ix)
		{
			FluidElement& phi = lab(ix,iy,iz);
			FluidElement& phiW = lab(ix-1,iy  ,iz  );
			FluidElement& phiE = lab(ix+1,iy  ,iz  );
			FluidElement& phiS = lab(ix  ,iy-1,iz  );
			FluidElement& phiN = lab(ix  ,iy+1,iz  );
			FluidElement& phiF = lab(ix  ,iy  ,iz-1);
			FluidElement& phiB = lab(ix  ,iy  ,iz+1);
			
			o(ix,iy,iz).tmpU = inv2h * (phiN.w-phiS.w) - inv2h * (phiB.v-phiF.v);
			o(ix,iy,iz).tmpV = inv2h * (phiB.u-phiB.u) - inv2h * (phiE.w-phiW.w);
			o(ix,iy,iz).tmpW = inv2h * (phiE.v-phiW.v) - inv2h * (phiN.u-phiS.u);
			o(ix,iy,iz).tmp = sqrt(o(ix,iy,iz).tmpU*o(ix,iy,iz).tmpU + o(ix,iy,iz).tmpV*o(ix,iy,iz).tmpV + o(ix,iy,iz).tmpW*o(ix,iy,iz).tmpW);
		}
	}
}

#endif /* OperatorVorticity_h */
