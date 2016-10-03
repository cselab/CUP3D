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
#include "InterpolationKernels.h"
#include "GenericOperator.h"
#include <cmath>

class OperatorAdvectionUpwind3rdOrder : public GenericLabOperator
{
private:
	const Real dt;
	const Real* const uInf;

public:
	OperatorAdvectionUpwind3rdOrder(const Real dt, const Real* const uInf)
: dt(dt), uInf(uInf)
	{
		stencil = StencilInfo(-2,-2,-2, 3,3,3, false, 3, 0,1,2);

		stencil_start[0] = -2;
		stencil_start[1] = -2;
		stencil_start[2] = -2;
		stencil_end[0] = 3;
		stencil_end[1] = 3;
		stencil_end[2] = 3;
	}

	~OperatorAdvectionUpwind3rdOrder() {}

	template <typename Lab, typename BlockType>
	void operator()(Lab & lab, const BlockInfo& info, BlockType& o) const
	{
#ifndef _RK2_
		const Real factor = -dt/(6.*info.h_gridpoint);
#else //perform half step
		const Real factor = -dt/(12.*info.h_gridpoint);
#endif

		for (int iz=0; iz<FluidBlock::sizeZ; ++iz)
			for (int iy=0; iy<FluidBlock::sizeY; ++iy)
				for (int ix=0; ix<FluidBlock::sizeX; ++ix)
				{
					const FluidElement& phi   = lab(ix  ,iy  ,iz  );
					const FluidElement& phiW  = lab(ix-1,iy  ,iz  );
					const FluidElement& phiE  = lab(ix+1,iy  ,iz  );
					const FluidElement& phiS  = lab(ix  ,iy-1,iz  );
					const FluidElement& phiN  = lab(ix  ,iy+1,iz  );
					const FluidElement& phiF  = lab(ix  ,iy  ,iz-1);
					const FluidElement& phiB  = lab(ix  ,iy  ,iz+1);
					const FluidElement& phiW2 = lab(ix-2,iy  ,iz  );
					const FluidElement& phiE2 = lab(ix+2,iy  ,iz  );
					const FluidElement& phiS2 = lab(ix  ,iy-2,iz  );
					const FluidElement& phiN2 = lab(ix  ,iy+2,iz  );
					const FluidElement& phiF2 = lab(ix  ,iy  ,iz-2);
					const FluidElement& phiB2 = lab(ix  ,iy  ,iz+2);

					const Real u3 = 3*phi.u;
					const Real v3 = 3*phi.v;
					const Real w3 = 3*phi.w;
					const Real u = phi.u + uInf[0];
					const Real v = phi.v + uInf[1];
					const Real w = phi.w + uInf[2];

					const Real dudx[2] = {  2*phiE.u + u3 - 6*phiW.u +   phiW2.u,
										   -  phiE2.u + 6*phiE.u - u3 - 2*phiW.u};
					const Real dvdx[2] = {  2*phiE.v + v3 - 6*phiW.v +   phiW2.v,
										   -  phiE2.v + 6*phiE.v - v3 - 2*phiW.v};
					const Real dwdx[2] = {  2*phiE.w + w3 - 6*phiW.w +   phiW2.w,
										   -  phiE2.w + 6*phiE.w - w3 - 2*phiW.w};

					const Real dudy[2] = {  2*phiN.u + u3 - 6*phiS.u +   phiS2.u,
										   -  phiN2.u + 6*phiN.u - u3 - 2*phiS.u};
					const Real dvdy[2] = {  2*phiN.v + v3 - 6*phiS.v +   phiS2.v,
										   -  phiN2.v + 6*phiN.v - v3 - 2*phiS.v};
					const Real dwdy[2] = {  2*phiN.w + w3 - 6*phiS.w +   phiS2.w,
										   -  phiN2.w + 6*phiN.w - w3 - 2*phiS.w};

					const Real dudz[2] = {  2*phiB.u + u3 - 6*phiF.u +   phiF2.u,
										   -  phiB2.u + 6*phiB.u - u3 - 2*phiF.u};
					const Real dvdz[2] = {  2*phiB.v + v3 - 6*phiF.v +   phiF2.v,
										   -  phiB2.v + 6*phiB.v - v3 - 2*phiF.v};
					const Real dwdz[2] = {  2*phiB.w + w3 - 6*phiF.w +   phiF2.w,
										   -  phiB2.w + 6*phiB.w - w3 - 2*phiF.w};

					const Real maxu = max(u,(Real)0);
					const Real maxv = max(v,(Real)0);
					const Real maxw = max(w,(Real)0);
					const Real minu = min(u,(Real)0);
					const Real minv = min(v,(Real)0);
					const Real minw = min(w,(Real)0);

					o(ix,iy,iz).tmpU = phi.u + factor*(maxu * dudx[0] + minu * dudx[1] +
													   maxv * dudy[0] + minv * dudy[1] +
													   maxw * dudz[0] + minw * dudz[1]);
					o(ix,iy,iz).tmpV = phi.v + factor*(maxu * dvdx[0] + minu * dvdx[1] +
													   maxv * dvdy[0] + minv * dvdy[1] +
													   maxw * dvdz[0] + minw * dvdz[1]);
					o(ix,iy,iz).tmpW = phi.w + factor*(maxu * dwdx[0] + minu * dwdx[1] +
													   maxv * dwdy[0] + minv * dwdy[1] +
													   maxw * dwdz[0] + minw * dwdz[1]);
				}
	}
};

#ifdef _RK2_
class OperatorAdvectionUpwind3rdOrderStage2 : public GenericLabOperator
{
private:
	const Real dt;
	const Real* const uInf;

public:
	OperatorAdvectionUpwind3rdOrderStage2(const Real dt, const Real* const uInf)
: dt(dt), uInf(uInf)
	{
		stencil = StencilInfo(-2,-2,-2, 3,3,3, false, 3, 5,6,7);
		stencil_start[0] = -2;
		stencil_start[1] = -2;
		stencil_start[2] = -2;
		stencil_end[0] = 3;
		stencil_end[1] = 3;
		stencil_end[2] = 3;
	}

	~OperatorAdvectionUpwind3rdOrderStage2() {}

	template <typename Lab, typename BlockType>
	void operator()(Lab & lab, const BlockInfo& info, BlockType& o) const
	{
		const Real factor = -dt/(6.*info.h_gridpoint);

		for (int iz=0; iz<FluidBlock::sizeZ; ++iz)
			for (int iy=0; iy<FluidBlock::sizeY; ++iy)
				for (int ix=0; ix<FluidBlock::sizeX; ++ix)
				{
					const FluidElement& phi  = lab(ix,iy,iz);
					const FluidElement& phiW = lab(ix-1,iy  ,iz  );
					const FluidElement& phiE = lab(ix+1,iy  ,iz  );
					const FluidElement& phiS = lab(ix  ,iy-1,iz  );
					const FluidElement& phiN = lab(ix  ,iy+1,iz  );
					const FluidElement& phiF = lab(ix  ,iy  ,iz-1);
					const FluidElement& phiB = lab(ix  ,iy  ,iz+1);
					const FluidElement& phiW2 = lab(ix-2,iy  ,iz  );
					const FluidElement& phiE2 = lab(ix+2,iy  ,iz  );
					const FluidElement& phiS2 = lab(ix  ,iy-2,iz  );
					const FluidElement& phiN2 = lab(ix  ,iy+2,iz  );
					const FluidElement& phiF2 = lab(ix  ,iy  ,iz-2);
					const FluidElement& phiB2 = lab(ix  ,iy  ,iz+2);
					const Real u3 = 3*phi.tmpU;
					const Real v3 = 3*phi.tmpV;
					const Real w3 = 3*phi.tmpW;

					const Real dudx[2] = {  2*phiE.tmpU + u3 - 6*phiW.tmpU +   phiW2.tmpU,
										   -  phiE2.tmpU + 6*phiE.tmpU - u3 - 2*phiW.tmpU};

					const Real dudy[2] = {  2*phiS.tmpU + u3 - 6*phiS.tmpU +   phiS2.tmpU,
										   -  phiN2.tmpU + 6*phiS.tmpU - u3 - 2*phiS.tmpU};

					const Real dudz[2] = {  2*phiB.tmpU + u3 - 6*phiF.tmpU +   phiF2.tmpU,
										   -  phiB2.tmpU + 6*phiB.tmpU - u3 - 2*phiF.tmpU};

					const Real dvdx[2] = {  2*phiE.tmpV + v3 - 6*phiW.tmpV +   phiW2.tmpV,
										   -  phiE2.tmpV + 6*phiE.tmpV - v3 - 2*phiW.tmpV};

					const Real dvdy[2] = {  2*phiS.tmpV + v3 - 6*phiS.tmpV +   phiS2.tmpV,
										   -  phiN2.tmpV + 6*phiS.tmpV - v3 - 2*phiS.tmpV};

					const Real dvdz[2] = {  2*phiB.tmpV + v3 - 6*phiF.tmpV +   phiF2.tmpV,
										   -  phiB2.tmpV + 6*phiB.tmpV - v3 - 2*phiF.tmpV};

					const Real dwdx[2] = {  2*phiE.tmpW + w3 - 6*phiW.tmpW +   phiW2.tmpW,
										   -  phiE2.tmpW + 6*phiE.tmpW - w3 - 2*phiW.tmpW};

					const Real dwdy[2] = {  2*phiS.tmpW + w3 - 6*phiS.tmpW +   phiS2.tmpW,
										   -  phiN2.tmpW + 6*phiS.tmpW - w3 - 2*phiS.tmpW};

					const Real dwdz[2] = {  2*phiB.tmpW + w3 - 6*phiF.tmpW +   phiF2.tmpW,
										   -  phiB2.tmpW + 6*phiB.tmpW - w3 - 2*phiF.tmpW};

					const Real u = phi.tmpU + uInf[0];
					const Real v = phi.tmpV + uInf[1];
					const Real w = phi.tmpW + uInf[2];
					const Real maxu = max(u,(Real)0);
					const Real maxv = max(v,(Real)0);
					const Real maxw = max(w,(Real)0);
					const Real minu = min(u,(Real)0);
					const Real minv = min(v,(Real)0);
					const Real minw = min(w,(Real)0);

					o(ix,iy,iz).u += factor*(maxu * dudx[0] + minu * dudx[1] +
											 maxv * dudy[0] + minv * dudy[1] +
											 maxw * dudz[0] + minw * dudz[1]);
					o(ix,iy,iz).v += factor*(maxu * dvdx[0] + minu * dvdx[1] +
										     maxv * dvdy[0] + minv * dvdy[1] +
										     maxw * dvdz[0] + minw * dvdz[1]);
					o(ix,iy,iz).w += factor*(maxu * dwdx[0] + minu * dwdx[1] +
											 maxv * dwdy[0] + minv * dwdy[1] +
											 maxw * dwdz[0] + minw * dwdz[1]);
				}
	}
};
#endif // _RK2_

template <typename Lab>
class CoordinatorAdvection : public GenericCoordinator
{
protected:
	const Real* const uInf;
	
	inline void update()
	{
		const int N = vInfo.size();
		
#pragma omp parallel for schedule(static)
		for(int i=0; i<N; i++) {
			BlockInfo info = vInfo[i];
			FluidBlock& b = *(FluidBlock*)info.ptrBlock;
			
			for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
				for(int iy=0; iy<FluidBlock::sizeY; ++iy)
					for(int ix=0; ix<FluidBlock::sizeX; ++ix)
					{
						b(ix,iy,iz).u = b(ix,iy,iz).tmpU;
						b(ix,iy,iz).v = b(ix,iy,iz).tmpV;
						b(ix,iy,iz).w = b(ix,iy,iz).tmpW;
					}
		}
	}
	
public:
	CoordinatorAdvection(const Real* const uInf, FluidGridMPI * grid)
	: GenericCoordinator(grid), uInf(uInf)
	{ }
	
	~CoordinatorAdvection()
	{ }
	
	void operator()(const Real dt)
	{
		check("advection - start");
		
		{
			OperatorAdvectionUpwind3rdOrder kernel(dt,uInf);
			compute(kernel);
		}
#ifdef _RK2_
		{
			OperatorAdvectionUpwind3rdOrderStage2 kernel(dt,uInf);
			compute(kernel);
		}
#else
		update();
#endif
		
		check("advection - end");
	}
	
	string getName()
	{
		return "Advection";
	}
};
#endif
