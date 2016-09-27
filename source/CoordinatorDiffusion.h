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
#include "GenericOperator.h"

class OperatorDiffusion : public GenericLabOperator
{
private:
	const Real mu;
	Real dt;

public:
	OperatorDiffusion(Real dt, Real mu) : mu(mu), dt(dt)
	{
		stencil = StencilInfo(-1,-1,-1, 2,2,2, false, 3, 0,1,2);
		stencil_start[0] = -1;
		stencil_start[1] = -1;
		stencil_start[2] = -1;
		stencil_end[0] = 2;
		stencil_end[1] = 2;
		stencil_end[2] = 2;
	}

	~OperatorDiffusion() {}

	template <typename Lab, typename BlockType>
	void operator()(Lab & lab, const BlockInfo& info, BlockType& o) const
	{
#ifdef _RK2_
		const Real fac = mu * 0.5 * dt / (info.h_gridpoint*info.h_gridpoint);
#else
		const Real fac = mu * dt / (info.h_gridpoint*info.h_gridpoint);
#endif // _RK2_

		for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
		for(int iy=0; iy<FluidBlock::sizeY; ++iy)
		for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
			FluidElement& phi  = lab(ix,iy,iz);
			FluidElement& phiW = lab(ix-1,iy  ,iz  );
			FluidElement& phiE = lab(ix+1,iy  ,iz  );
			FluidElement& phiS = lab(ix  ,iy-1,iz  );
			FluidElement& phiN = lab(ix  ,iy+1,iz  );
			FluidElement& phiF = lab(ix  ,iy  ,iz-1);
			FluidElement& phiB = lab(ix  ,iy  ,iz+1);

			o(ix,iy,iz).tmpU = phi.u + fac * (phiN.u + phiS.u + phiE.u + phiW.u + phiF.u + phiB.u - phi.u*6.);
			o(ix,iy,iz).tmpV = phi.v + fac * (phiN.v + phiS.v + phiE.v + phiW.v + phiF.v + phiB.v - phi.v*6.);
			o(ix,iy,iz).tmpW = phi.w + fac * (phiN.w + phiS.w + phiE.w + phiW.w + phiF.w + phiB.w - phi.w*6.);
		}
	}
};

#ifdef _RK2_
class OperatorDiffusionStage2 : public GenericLabOperator
{
private:
	const Real mu;
	Real dt;

public:
	OperatorDiffusionStage2(Real dt, Real mu) : mu(mu), dt(dt)
	{
		stencil = StencilInfo(-1,-1,-1, 2,2,2, false, 3, 5,6,7);

		stencil_start[0] = -1;
		stencil_start[1] = -1;
		stencil_start[2] = -1;
		stencil_end[0] = 2;
		stencil_end[1] = 2;
		stencil_end[2] = 2;
	}

	~OperatorDiffusionStage2() {}

	template <typename Lab, typename BlockType>
	void operator()(Lab & lab, const BlockInfo& info, BlockType& o) const
	{
		const Real fac = mu * dt / (info.h_gridpoint*info.h_gridpoint);

		for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
		for(int iy=0; iy<FluidBlock::sizeY; ++iy)
		for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
			FluidElement& phi = lab(ix,iy,iz);
			FluidElement& phiW = lab(ix-1,iy  ,iz  );
			FluidElement& phiE = lab(ix+1,iy  ,iz  );
			FluidElement& phiS = lab(ix  ,iy-1,iz  );
			FluidElement& phiN = lab(ix  ,iy+1,iz  );
			FluidElement& phiF = lab(ix  ,iy  ,iz-1);
			FluidElement& phiB = lab(ix  ,iy  ,iz+1);

			o(ix,iy,iz).u += fac * (phiN.tmpU + phiS.tmpU + phiE.tmpU + phiW.tmpU + phiF.tmpU + phiB.tmpU - phi.tmpU*6.);
			o(ix,iy,iz).v += fac * (phiN.tmpV + phiS.tmpV + phiE.tmpV + phiW.tmpV + phiF.tmpV + phiB.tmpV - phi.tmpV*6.);
			o(ix,iy,iz).w += fac * (phiN.tmpW + phiS.tmpW + phiE.tmpW + phiW.tmpW + phiF.tmpW + phiB.tmpW - phi.tmpW*6.);
		}
	}
};
#endif // _RK2_

template <typename Lab>
class CoordinatorDiffusion : public GenericCoordinator
{
protected:
	const Real coeff;
	
	inline void update()
	{
#pragma omp parallel for schedule(static)
		for(int i=0; i<vInfo.size(); i++) {
			BlockInfo info = vInfo[i];
			FluidBlock& b = *(FluidBlock*)info.ptrBlock;
			
			for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
			for(int iy=0; iy<FluidBlock::sizeY; ++iy)
			for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
				b(ix,iy,iz).u = b(ix,iy,iz).tmpU;
				b(ix,iy,iz).v = b(ix,iy,iz).tmpV;
				b(ix,iy,iz).w = b(ix,iy,iz).tmpW;
			}
		}
	}
	
	inline void diffuse(const Real dt, const int stage)
	{

	}
	
public:
	CoordinatorDiffusion(const Real coeff, FluidGridMPI * grid) : GenericCoordinator(grid), coeff(coeff)
	{
	}
	
	void operator()(const Real dt)
	{
		check("diffusion - start");
		
		{
			OperatorDiffusion kernel(dt, coeff);
			compute(kernel);
		}
		
#ifdef _RK2_
		{
			OperatorDiffusionStage2 kernel(dt, coeff);
			compute(kernel);
		}
#else
		update();
#endif
		
		check("diffusion - end");
	}
	
	string getName()
	{
		return "Diffusion";
	}
};

#endif
