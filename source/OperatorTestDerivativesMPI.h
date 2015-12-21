//
//  OperatorTestDerivativesMPI.h
//  CubismUP_3D
//
//  Created by Christian Conti on 12/18/15.
//  Copyright © 2015 ETHZ. All rights reserved.
//

#ifndef OperatorTestDerivativesMPI_h
#define OperatorTestDerivativesMPI_h

class OperatorTestDerivativesMPI : public GenericLabOperator
{
private:
	double dt;
	
public:
	OperatorTestDerivativesMPI(double dt) : dt(dt)
	{
		stencil = StencilInfo(-1,-1,-1, 2,2,2, false, 3, 1,2,3);
		
		stencil_start[0] = -1;
		stencil_start[1] = -1;
		stencil_start[2] = -1;
		stencil_end[0] = 2;
		stencil_end[1] = 2;
		stencil_end[2] = 2;
	}
	
	~OperatorTestDerivativesMPI() {}
	
	template <typename Lab, typename BlockType>
	void operator()(Lab & lab, const BlockInfo& info, BlockType& o) const
	{
		const double prefactor = 1. / (info.h_gridpoint*info.h_gridpoint);
		
		for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
			for(int iy=0; iy<FluidBlock::sizeY; ++iy)
				for(int ix=0; ix<FluidBlock::sizeX; ++ix)
				{
					FluidElement& phi = lab(ix,iy,iz);
					FluidElement& phiE = lab(ix+1,iy  ,iz  );
					FluidElement& phiW = lab(ix-1,iy  ,iz  );
					FluidElement& phiN = lab(ix  ,iy+1,iz  );
					FluidElement& phiS = lab(ix  ,iy-1,iz  );
					FluidElement& phiB = lab(ix  ,iy  ,iz+1);
					FluidElement& phiF = lab(ix  ,iy  ,iz-1);
					
					o(ix,iy,iz).chi = prefactor * (phiN.u + phiS.u + phiE.u + phiW.u + phiF.u + phiB.u - phi.u*6.);
					o(ix,iy,iz).p   = prefactor * (phiN.v + phiS.v + phiE.v + phiW.v + phiF.v + phiB.v - phi.v*6.);
					o(ix,iy,iz).tmp = prefactor * (phiN.w + phiS.w + phiE.w + phiW.w + phiF.w + phiB.w - phi.w*6.);
				}
	}
};

#endif /* OperatorTestDerivativesMPI_h */
