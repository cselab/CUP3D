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

#ifdef _ISPC_
#include "Advection_ispc.h"

using namespace ispc;

class OperatorAdvectionUpwind3rdOrderISPC : public GenericLabOperator
{
private:
	double dt;
	const int stage;
	Real *uBody, *vBody, *wBody;
	
	Real *uLab, *vLab, *wLab;
	Real *u, *v, *w;
	Real *uTmp, *vTmp, *wTmp;
	
	enum {
		_BSX2_ = _BSX_+4,
		_BSY2_ = _BSY_+4,
		_BSZ2_ = _BSZ_+4,
		SLICESIZE = _BSX_*_BSY_,
		SLICESIZE2 = (_BSX_+4)*(_BSY_+4),
		SIZE = _BSX_*_BSY_*_BSZ_,
		SIZE2 = (_BSX_+4)*(_BSY_+4)*(_BSZ_+4)
	};
	
public:
	OperatorAdvectionUpwind3rdOrderISPC(double dt, Real * uBody, Real * vBody, Real * wBody, Real * u, Real * v, Real * w, Real * uTmp, Real * vTmp, Real * wTmp, Real * uLab, Real * vLab, Real * wLab, const int stage) : dt(dt), uBody(uBody), vBody(vBody), wBody(wBody), stage(stage), u(u), v(v), w(w), uTmp(uTmp), vTmp(vTmp), wTmp(wTmp), uLab(uLab), vLab(vLab), wLab(wLab)
	{
		stencil = StencilInfo(-2,-2,-2, 3,3,3, false, 8, 0,1,2,3,7,8,9,10);
		
		stencil_start[0] = -2;
		stencil_start[1] = -2;
		stencil_start[2] = -2;
		
		stencil_end[0] = 3;
		stencil_end[1] = 3;
		stencil_end[2] = 3;
	}
	
	OperatorAdvectionUpwind3rdOrderISPC(double dt, Real * u, Real * v, Real * w, Real * uTmp, Real * vTmp, Real * wTmp, Real * uLab, Real * vLab, Real * wLab, const int stage) : dt(dt), uBody(NULL), vBody(NULL), wBody(NULL), stage(stage), u(u), v(v), w(w), uTmp(uTmp), vTmp(vTmp), wTmp(wTmp), uLab(uLab), vLab(vLab), wLab(wLab)
	{
		stencil = StencilInfo(-2,-2,-2, 3,3,3, false, 8, 0,1,2,3,7,8,9,10);
		
		stencil_start[0] = -2;
		stencil_start[1] = -2;
		stencil_start[2] = -2;
		
		stencil_end[0] = 3;
		stencil_end[1] = 3;
		stencil_end[2] = 3;
	}
	
	~OperatorAdvectionUpwind3rdOrderISPC()
	{
	}
	
	template <typename Lab, typename BlockType>
	void operator()(Lab & lab, const BlockInfo& info, BlockType& o) const
	{
		const Real factor = -dt/(6.*info.h_gridpoint);
		
		const int tid = omp_get_thread_num();
		/*
		 // AoS -> SoA
		for (int iz=0; iz<_BSZ2_; ++iz)
			for (int iy=0; iy<_BSY2_; ++iy)
				for (int ix=0; ix<_BSX2_; ++ix)
				{
					const int idx = ix + iy*_BSX2_ + iz*SLICESIZE2;
					assert(idx>=0 && idx<SIZE2);
					
					uLab[tid*SIZE2 + idx] = lab(ix-2,iy-2,iz-2).u;
					vLab[tid*SIZE2 + idx] = lab(ix-2,iy-2,iz-2).v;
					wLab[tid*SIZE2 + idx] = lab(ix-2,iy-2,iz-2).w;
				}
		 
		for (int iz=0; iz<_BSZ_; ++iz)
			for (int iy=0; iy<_BSY_; ++iy)
				for (int ix=0; ix<_BSX_; ++ix)
				{
					const int idx = ix + iy*_BSX_ + iz*SLICESIZE;
					assert(idx>=0 && idx<SIZE);
					
					u[tid*SIZE + idx] = o(ix,iy,iz).u;
					v[tid*SIZE + idx] = o(ix,iy,iz).v;
					w[tid*SIZE + idx] = o(ix,iy,iz).w;
				}
		*/
		// ISPC
		const int offset = tid*SIZE2 + 2*(1 + _BSX2_ + SLICESIZE2);
#ifdef _MOVING_FRAME_
		advect_ispc(u,v,w,uTmp,vTmp,wTmp,uLab+offset,vLab+offset,wLab+offset,*uBody,*vBody,*wBody,factor);
#else
		advect_ispc(u,v,w,uTmp,vTmp,wTmp,uLab+offset,vLab+offset,wLab+offset,factor);
#endif
		/*
		// SoA -> AoS
		for (int iz=0; iz<_BSZ_; ++iz)
			for (int iy=0; iy<_BSY_; ++iy)
				for (int ix=0; ix<_BSX_; ++ix)
				{
					const int idx = ix + iy*_BSX_ + iz*SLICESIZE;
					assert(idx>=0 && idx<SIZE);
					
					o(ix,iy,iz).tmpU = uTmp[tid*SIZE + idx];
					o(ix,iy,iz).tmpV = vTmp[tid*SIZE + idx];
					o(ix,iy,iz).tmpW = wTmp[tid*SIZE + idx];
				}
		 */
	}
	
};
#endif

class OperatorAdvectionUpwind3rdOrder : public GenericLabOperator
{
private:
	double dt;
	const int stage;
	Real *uBody, *vBody, *wBody;
	
public:
	OperatorAdvectionUpwind3rdOrder(double dt, Real * uBody, Real * vBody, Real * wBody, const int stage) : dt(dt), uBody(uBody), vBody(vBody), wBody(wBody), stage(stage)
	{
		stencil = StencilInfo(-2,-2,-2, 3,3,3, false, 8, 0,1,2,3,7,8,9,10);
		
		stencil_start[0] = -2;
		stencil_start[1] = -2;
		stencil_start[2] = -2;
		
		stencil_end[0] = 3;
		stencil_end[1] = 3;
		stencil_end[2] = 3;
	}
	
	OperatorAdvectionUpwind3rdOrder(double dt, const int stage) : dt(dt), uBody(NULL), vBody(NULL), wBody(NULL), stage(stage)
	{
		stencil = StencilInfo(-2,-2,-2, 3,3,3, false, 8, 0,1,2,3,7,8,9,10);
		
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
#else
		const Real factor = -dt*(stage==0 ? .5 : 1)/(6.*info.h_gridpoint);
#endif
		
		if (stage==0)
			for (int iz=0; iz<FluidBlock::sizeZ; ++iz)
				for (int iy=0; iy<FluidBlock::sizeY; ++iy)
					for (int ix=0; ix<FluidBlock::sizeX; ++ix)
					{
						FluidElement& phi  = lab(ix,iy,iz);
						FluidElement& phiW = lab(ix-1,iy  ,iz  );
						FluidElement& phiE = lab(ix+1,iy  ,iz  );
						FluidElement& phiS = lab(ix  ,iy-1,iz  );
						FluidElement& phiN = lab(ix  ,iy+1,iz  );
						FluidElement& phiF = lab(ix  ,iy  ,iz-1);
						FluidElement& phiB = lab(ix  ,iy  ,iz+1);
						FluidElement& phiW2 = lab(ix-2,iy  ,iz  );
						FluidElement& phiE2 = lab(ix+2,iy  ,iz  );
						FluidElement& phiS2 = lab(ix  ,iy-2,iz  );
						FluidElement& phiN2 = lab(ix  ,iy+2,iz  );
						FluidElement& phiF2 = lab(ix  ,iy  ,iz-2);
						FluidElement& phiB2 = lab(ix  ,iy  ,iz+2);
						
						const Real u3 = 3*phi.u;
						const Real v3 = 3*phi.v;
						const Real w3 = 3*phi.w;
						
						const Real dudx[2] = {  2*phiE.u + u3 - 6*phiW.u +   phiW2.u,
											   -  phiE2.u + 6*phiE.u - u3 - 2*phiW.u};
						
						const Real dudy[2] = {  2*phiN.u + u3 - 6*phiS.u +   phiS2.u,
											   -  phiN2.u + 6*phiN.u - u3 - 2*phiS.u};
						
						const Real dudz[2] = {  2*phiB.u + u3 - 6*phiF.u +   phiF2.u,
											   -  phiB2.u + 6*phiB.u - u3 - 2*phiF.u};
						
						const Real dvdx[2] = {  2*phiE.v + v3 - 6*phiW.v +   phiW2.v,
											   -  phiE2.v + 6*phiE.v - v3 - 2*phiW.v};
						
						const Real dvdy[2] = {  2*phiN.v + v3 - 6*phiS.v +   phiS2.v,
											   -  phiN2.v + 6*phiN.v - v3 - 2*phiS.v};
						
						const Real dvdz[2] = {  2*phiB.v + v3 - 6*phiF.v +   phiF2.v,
											   -  phiB2.v + 6*phiB.v - v3 - 2*phiF.v};
						
						const Real dwdx[2] = {  2*phiE.w + w3 - 6*phiW.w +   phiW2.w,
											   -  phiE2.w + 6*phiE.w - w3 - 2*phiW.w};
						
						const Real dwdy[2] = {  2*phiN.w + w3 - 6*phiS.w +   phiS2.w,
											   -  phiN2.w + 6*phiN.w - w3 - 2*phiS.w};
						
						const Real dwdz[2] = {  2*phiB.w + w3 - 6*phiF.w +   phiF2.w,
											   -  phiB2.w + 6*phiB.w - w3 - 2*phiF.w};
						
#ifndef _MOVING_FRAME_
						const Real u = o(ix,iy,iz).u;
						const Real v = o(ix,iy,iz).v;
						const Real w = o(ix,iy,iz).w;
#else
						const Real u = o(ix,iy,iz).u - *uBody;
						const Real v = o(ix,iy,iz).v - *vBody;
						const Real w = o(ix,iy,iz).w - *wBody;
#endif
						
						const Real maxu = max(u,(Real)0);
						const Real maxv = max(v,(Real)0);
						const Real maxw = max(w,(Real)0);
						const Real minu = min(u,(Real)0);
						const Real minv = min(v,(Real)0);
						const Real minw = min(w,(Real)0);
						
						o(ix,iy,iz).tmpU = o(ix,iy,iz).u + factor*(maxu * dudx[0] + minu * dudx[1] +
																   maxv * dudy[0] + minv * dudy[1] +
																   maxw * dudz[0] + minw * dudz[1]);
						o(ix,iy,iz).tmpV = o(ix,iy,iz).v + factor*(maxu * dvdx[0] + minu * dvdx[1] +
																   maxv * dvdy[0] + minv * dvdy[1] +
																   maxw * dvdz[0] + minw * dvdz[1]);
						o(ix,iy,iz).tmpW = o(ix,iy,iz).w + factor*(maxu * dwdx[0] + minu * dwdx[1] +
																   maxv * dwdy[0] + minv * dwdy[1] +
																   maxw * dwdz[0] + minw * dwdz[1]);
#ifndef _RK2_
#ifdef _MULTIPHASE_
						const Real r3 = 3*phi.rho;
						
						const Real drdx[2] = {  2*phiE.rho + r3 - 6*phiW.rho +   phiW2.rho,
											   -  phiE2.rho + 6*phiE.rho - r3 - 2*phiW.rho};
						
						const Real drdy[2] = {  2*phiN.rho + r3 - 6*phiS.rho +   phiS2.rho,
											   -  phiN2.rho + 6*phiN.rho - r3 - 2*phiS.rho};
						
						const Real drdz[2] = {  2*phiB.rho + r3 - 6*phiF.rho +   phiF2.rho,
											   -  phiB2.rho + 6*phiB.rho - r3 - 2*phiF.rho};
						
						const Real r = o(ix,iy,iz).rho;
						
						o(ix,iy,iz).tmp  = r + factor*(maxu * drdx[0] + minu * drdx[1] +
													   maxv * drdy[0] + minv * drdy[1] +
													   maxw * drdz[0] + minw * drdz[1]);
#endif // _MULTIPHASE_
#endif // _RK2_
					}
#ifdef _RK2_
		else if (stage==1)
			for (int iz=0; iz<FluidBlock::sizeZ; ++iz)
				for (int iy=0; iy<FluidBlock::sizeY; ++iy)
					for (int ix=0; ix<FluidBlock::sizeX; ++ix)
					{
						FluidElement& phi  = lab(ix,iy,iz);
						FluidElement& phiW = lab(ix-1,iy  ,iz  );
						FluidElement& phiE = lab(ix+1,iy  ,iz  );
						FluidElement& phiS = lab(ix  ,iy-1,iz  );
						FluidElement& phiN = lab(ix  ,iy+1,iz  );
						FluidElement& phiF = lab(ix  ,iy  ,iz-1);
						FluidElement& phiB = lab(ix  ,iy  ,iz+1);
						FluidElement& phiW2 = lab(ix-2,iy  ,iz  );
						FluidElement& phiE2 = lab(ix+2,iy  ,iz  );
						FluidElement& phiS2 = lab(ix  ,iy-2,iz  );
						FluidElement& phiN2 = lab(ix  ,iy+2,iz  );
						FluidElement& phiF2 = lab(ix  ,iy  ,iz-2);
						FluidElement& phiB2 = lab(ix  ,iy  ,iz+2);
						
						const Real u3 = 3*phi.tmpU;
						const Real v3 = 3*phi.tmpV;
						const Real w3 = 3*phi.tmpW;
						const Real r3 = 3*phi.tmp;
						
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
						
#ifndef _MOVING_FRAME_
						const Real u = phi.tmpU;
						const Real v = phi.tmpV;
						const Real w = phi.tmpW;
#else
						const Real u = phi.tmpU - *uBody;
						const Real v = phi.tmpV - *vBody;
						const Real w = phi.tmpW - *wBody;
#endif
						
						const Real maxu = max(u,(Real)0);
						const Real maxv = max(v,(Real)0);
						const Real maxw = max(w,(Real)0);
						const Real minu = min(u,(Real)0);
						const Real minv = min(v,(Real)0);
						const Real minw = min(w,(Real)0);
						
						o(ix,iy,iz).tmpU = o(ix,iy,iz).u + factor*(maxu * dudx[0] + minu * dudx[1] +
																   maxv * dudy[0] + minv * dudy[1] +
																   maxw * dudz[0] + minw * dudz[1]);
						o(ix,iy,iz).tmpV = o(ix,iy,iz).v + factor*(maxu * dvdx[0] + minu * dvdx[1] +
																   maxv * dvdy[0] + minv * dvdy[1] +
																   maxw * dvdz[0] + minw * dvdz[1]);
						o(ix,iy,iz).tmpW = o(ix,iy,iz).w + factor*(maxu * dwdx[0] + minu * dwdx[1] +
																   maxv * dwdy[0] + minv * dwdy[1] +
																   maxw * dwdz[0] + minw * dwdz[1]);
#ifdef _MULTIPHASE_
						const Real drdx[2] = {  2*phiE.rho + 3*phi.rho - 6*phiW.rho +   phiW2.rho,
											   -  phiE2.rho + 6*phiE.rho - 3*phi.rho - 2*phiW.rho};
						
						const Real drdy[2] = {  2*phiS.rho + 3*phi.rho - 6*phiS.rho +   phiS2.rho,
											   -  phiN2.rho + 6*phiS.rho - 3*phi.rho - 2*phiS.rho};
						
						const Real drdz[2] = {  2*phiB.rho + 3*phi.rho - 6*phiF.rho +   phiF2.rho,
											   -  phiB2.rho + 6*phiB.rho - 3*phi.rho - 2*phiF.rho};
						
						o(ix,iy,iz).tmp = phi.rho + factor*(maxu * drdx[0] + minu * drdx[1] +
															maxv * drdy[0] + minv * drdy[1] +
															maxw * drdz[0] + minw * drdz[1]);
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
		stencil = StencilInfo(-2,-2,-2, 3,3,3, false, 8, 0,1,2,3,7,8,9,10);
		
		stencil_start[0] = -2;
		stencil_start[1] = -2;
		stencil_start[2] = -2;
		
		stencil_end[0] = 3;
		stencil_end[1] = 3;
		stencil_end[2] = 3;
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
			for (int iz=0; iz<FluidBlock::sizeZ; ++iz)
				for (int iy=0; iy<FluidBlock::sizeY; ++iy)
					for (int ix=0; ix<FluidBlock::sizeX; ++ix)
					{
						double p[3];
						info.pos(p, ix, iy);
						
						const Real u = o(ix,iy).u;
						const Real v = o(ix,iy).v;
						const Real w = o(ix,iy).w;
						const Real r = o(ix,iy).rho;
#ifndef _RK2_
						const Real drdx[2] = {  2*lab(ix+1,iy  ,iz  ).rho + 3*lab(ix  ,iy  ,iz  ).rho - 6*lab(ix-1,iy  ,iz  ).rho +   lab(ix-2,iy  ,iz  ).rho,
							-  lab(ix+2,iy  ,iz  ).rho + 6*lab(ix+1,iy  ,iz  ).rho - 3*lab(ix  ,iy  ,iz  ).rho - 2*lab(ix-1,iy  ,iz  ).rho};
						
						const Real drdy[2] = {  2*lab(ix  ,iy+1,iz  ).rho + 3*lab(ix  ,iy  ,iz  ).rho - 6*lab(ix  ,iy-1,iz  ).rho +   lab(ix  ,iy-2,iz  ).rho,
							-  lab(ix  ,iy+2,iz  ).rho + 6*lab(ix  ,iy+1,iz  ).rho - 3*lab(ix  ,iy  ,iz  ).rho - 2*lab(ix  ,iy-1,iz  ).rho};
						
						const Real drdz[2] = {  2*lab(ix  ,iy  ,iz+1).rho + 3*lab(ix  ,iy  ,iz  ).rho - 6*lab(ix  ,iy  ,iz-1).rho +   lab(ix  ,iy  ,iz-2).rho,
							-  lab(ix  ,iy  ,iz+2).rho + 6*lab(ix  ,iy  ,iz+1).rho - 3*lab(ix  ,iy  ,iz  ).rho - 2*lab(ix  ,iy  ,iz-1).rho};
						
						o(ix,iy,iz).tmp  = r + factor*(max(u,(Real)0) * drdx[0] + min(u,(Real)0) * drdx[1] +
													   max(v,(Real)0) * drdy[0] + min(v,(Real)0) * drdy[1] +
													   max(w,(Real)0) * drdz[0] + min(w,(Real)0) * drdz[1]);
#else
						const Real dudx[2] = {  2*lab(ix+1,iy  ,iz  ).u + 3*lab(ix  ,iy  ,iz  ).u - 6*lab(ix-1,iy  ,iz  ).u +   lab(ix-2,iy  ,iz  ).u,
							-  lab(ix+2,iy  ,iz  ).u + 6*lab(ix+1,iy  ,iz  ).u - 3*lab(ix  ,iy  ,iz  ).u - 2*lab(ix-1,iy  ,iz  ).u};
						
						const Real dudy[2] = {  2*lab(ix  ,iy+1,iz  ).u + 3*lab(ix  ,iy  ,iz  ).u - 6*lab(ix  ,iy-1,iz  ).u +   lab(ix  ,iy-2,iz  ).u,
							-  lab(ix  ,iy+2,iz  ).u + 6*lab(ix  ,iy+1,iz  ).u - 3*lab(ix  ,iy  ,iz  ).u - 2*lab(ix  ,iy-1,iz  ).u};
						
						const Real dudz[2] = {  2*lab(ix  ,iy  ,iz+1).u + 3*lab(ix  ,iy  ,iz  ).u - 6*lab(ix  ,iy  ,iz-1).u +   lab(ix  ,iy  ,iz-2).u,
							-  lab(ix  ,iy  ,iz+2).u + 6*lab(ix  ,iy  ,iz+1).u - 3*lab(ix  ,iy  ,iz  ).u - 2*lab(ix  ,iy  ,iz-1).u};
						
						const Real dvdx[2] = {  2*lab(ix+1,iy  ,iz  ).v + 3*lab(ix  ,iy  ,iz  ).v - 6*lab(ix-1,iy  ,iz  ).v +   lab(ix-2,iy  ,iz  ).v,
							-  lab(ix+2,iy  ,iz  ).v + 6*lab(ix+1,iy  ,iz  ).v - 3*lab(ix  ,iy  ,iz  ).v - 2*lab(ix-1,iy  ,iz  ).v};
						
						const Real dvdy[2] = {  2*lab(ix  ,iy+1,iz  ).v + 3*lab(ix  ,iy  ,iz  ).v - 6*lab(ix  ,iy-1,iz  ).v +   lab(ix  ,iy-2,iz  ).v,
							-  lab(ix  ,iy+2,iz  ).v + 6*lab(ix  ,iy+1,iz  ).v - 3*lab(ix  ,iy  ,iz  ).v - 2*lab(ix  ,iy-1,iz  ).v};
						
						const Real dvdz[2] = {  2*lab(ix  ,iy  ,iz+1).v + 3*lab(ix  ,iy  ,iz  ).v - 6*lab(ix  ,iy  ,iz-1).v +   lab(ix  ,iy  ,iz-2).v,
							-  lab(ix  ,iy  ,iz+2).v + 6*lab(ix  ,iy  ,iz+1).v - 3*lab(ix  ,iy  ,iz  ).v - 2*lab(ix  ,iy  ,iz-1).v};
						
						const Real dwdx[2] = {  2*lab(ix+1,iy  ,iz  ).w + 3*lab(ix  ,iy  ,iz  ).w - 6*lab(ix-1,iy  ,iz  ).w +   lab(ix-2,iy  ,iz  ).w,
							-  lab(ix+2,iy  ,iz  ).w + 6*lab(ix+1,iy  ,iz  ).w - 3*lab(ix  ,iy  ,iz  ).w - 2*lab(ix-1,iy  ,iz  ).w};
						
						const Real dwdy[2] = {  2*lab(ix  ,iy+1,iz  ).w + 3*lab(ix  ,iy  ,iz  ).w - 6*lab(ix  ,iy-1,iz  ).w +   lab(ix  ,iy-2,iz  ).w,
							-  lab(ix  ,iy+2,iz  ).w + 6*lab(ix  ,iy+1,iz  ).w - 3*lab(ix  ,iy  ,iz  ).w - 2*lab(ix  ,iy-1,iz  ).w};
						
						const Real dwdz[2] = {  2*lab(ix  ,iy  ,iz+1).w + 3*lab(ix  ,iy  ,iz  ).w - 6*lab(ix  ,iy  ,iz-1).w +   lab(ix  ,iy  ,iz-2).w,
							-  lab(ix  ,iy  ,iz+2).w + 6*lab(ix  ,iy  ,iz+1).w - 3*lab(ix  ,iy  ,iz  ).w - 2*lab(ix  ,iy  ,iz-1).w};
						
						o(ix,iy,iz).tmpU = u + factor*(max(u,(Real)0) * dudx[0] + min(u,(Real)0) * dudx[1] +
													   max(v,(Real)0) * dudy[0] + min(v,(Real)0) * dudy[1] +
													   max(w,(Real)0) * dudz[0] + min(w,(Real)0) * dudz[1]);
						o(ix,iy,iz).tmpV = v + factor*(max(u,(Real)0) * dvdx[0] + min(u,(Real)0) * dvdx[1] +
													   max(v,(Real)0) * dvdy[0] + min(v,(Real)0) * dvdy[1] +
													   max(w,(Real)0) * dvdz[0] + min(w,(Real)0) * dvdz[1]);
						o(ix,iy,iz).tmpW = w + factor*(max(u,(Real)0) * dwdx[0] + min(u,(Real)0) * dwdx[1] +
													   max(v,(Real)0) * dwdy[0] + min(v,(Real)0) * dwdy[1] +
													   max(w,(Real)0) * dwdz[0] + min(w,(Real)0) * dwdz[1]);
#endif
					}
#ifdef _RK2_
		else if (stage==1)
			for (int iz=0; iz<FluidBlock::sizeZ; ++iz)
				for (int iy=0; iy<FluidBlock::sizeY; ++iy)
					for (int ix=0; ix<FluidBlock::sizeX; ++ix)
					{
						double p[3];
						info.pos(p, ix, iy);
						
						const Real drdx[2] = {  2*lab(ix+1,iy  ,iz  ).tmp + 3*lab(ix  ,iy  ,iz  ).tmp - 6*lab(ix-1,iy  ,iz  ).tmp +   lab(ix-2,iy  ,iz  ).tmp,
							-  lab(ix+2,iy  ,iz  ).tmp + 6*lab(ix+1,iy  ,iz  ).tmp - 3*lab(ix  ,iy  ,iz  ).tmp - 2*lab(ix-1,iy  ,iz  ).tmp};
						
						const Real drdy[2] = {  2*lab(ix  ,iy+1,iz  ).tmp + 3*lab(ix  ,iy  ,iz  ).tmp - 6*lab(ix  ,iy-1,iz  ).tmp +   lab(ix  ,iy-2,iz  ).tmp,
							-  lab(ix  ,iy+2,iz  ).tmp + 6*lab(ix  ,iy+1,iz  ).tmp - 3*lab(ix  ,iy  ,iz  ).tmp - 2*lab(ix  ,iy-1,iz  ).tmp};
						
						const Real drdz[2] = {  2*lab(ix  ,iy  ,iz+1).tmp + 3*lab(ix  ,iy  ,iz  ).tmp - 6*lab(ix  ,iy  ,iz-1).tmp +   lab(ix  ,iy  ,iz-2).tmp,
							-  lab(ix  ,iy  ,iz+2).tmp + 6*lab(ix  ,iy  ,iz+1).tmp - 3*lab(ix  ,iy  ,iz  ).tmp - 2*lab(ix  ,iy  ,iz-1).tmp};
						
						const Real u = lab(ix,iy).tmpU;
						const Real v = lab(ix,iy).tmpV;
						const Real w = lab(ix,iy).tmpW;
						const Real r = lab(ix,iy).rho;
						
						o(ix,iy,iz).tmp = r + factor*(max(u,(Real)0) * drdx[0] + min(u,(Real)0) * drdx[1] +
													  max(v,(Real)0) * drdy[0] + min(v,(Real)0) * drdy[1] +
													  max(w,(Real)0) * drdz[0] + min(w,(Real)0) * drdz[1]);
					}
#endif
	}
};

#endif
