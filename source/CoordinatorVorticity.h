//
//  CoordinatorVorticity.h
//  CubismUP_3D
//
//  Created by Christian Conti on 12/8/15.
//  Copyright © 2015 ETHZ. All rights reserved.
//

#ifndef CubismUP_3D_CoordinatorVorticity_h
#define CubismUP_3D_CoordinatorVorticity_h

#include "GenericCoordinator.h"
#include "GenericOperator.h"

class OperatorVorticity : public GenericLabOperator
{
private:

public:
	OperatorVorticity()
	{
		stencil = StencilInfo(-1,-1,-1, 2,2,2, false, 3, 1,2,3);

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
				for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
					FluidElement& phiW = lab(ix-1,iy  ,iz  );
					FluidElement& phiE = lab(ix+1,iy  ,iz  );
					FluidElement& phiS = lab(ix  ,iy-1,iz  );
					FluidElement& phiN = lab(ix  ,iy+1,iz  );
					FluidElement& phiF = lab(ix  ,iy  ,iz-1);
					FluidElement& phiB = lab(ix  ,iy  ,iz+1);

					o(ix,iy,iz).tmpU = inv2h * (phiN.w-phiS.w) - inv2h * (phiB.v-phiF.v);
					o(ix,iy,iz).tmpV = inv2h * (phiB.u-phiF.u) - inv2h * (phiE.w-phiW.w);
					o(ix,iy,iz).tmpW = inv2h * (phiE.v-phiW.v) - inv2h * (phiN.u-phiS.u);
				}
	}
};

template <typename Lab>
class CoordinatorVorticity : public GenericCoordinator
{
protected:
	
	inline void reset()
	{
		const int N = vInfo.size();
		
#pragma omp parallel for schedule(static)
		for(int i=0; i<N; i++) {
			BlockInfo info = vInfo[i];
			FluidBlock& b = *(FluidBlock*)info.ptrBlock;
			
			for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
				for(int iy=0; iy<FluidBlock::sizeY; ++iy)
					for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
						b(ix,iy,iz).tmpU = 0;
						b(ix,iy,iz).tmpV = 0;
						b(ix,iy,iz).tmpW = 0;
					}
		}
	};
	
public:
	CoordinatorVorticity(FluidGridMPI * grid) : GenericCoordinator(grid)
	{
	}
	
	void operator()(const double dt)
	{
		check("vorticity - start");
		
		//reset();

		OperatorVorticity kernel;
		compute(kernel);
		
		check("vorticity - end");
	}
	
	string getName()
	{
		return "Vorticity";
	}
};


#endif /* CoordinatorVorticity_h */
