//
//  OperatorDivergence.h
//  CubismUP_3D
//
//	Operates on
//		divU
//
//  Created by Christian Conti on 1/9/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef CubismUP_3D_OperatorDivergence_h
#define CubismUP_3D_OperatorDivergence_h

#include "Layer.h"
#include "GenericOperator.h"

class OperatorDivergenceLayer : public GenericLabOperator
{
private:
	Layer & divergence;
	
public:
	OperatorDivergenceLayer(Layer & divergence) : divergence(divergence)
	{
		stencil_start[0] = -1;
		stencil_start[1] = -1;
		stencil_start[2] = 0;
		stencil_end[0] = 2;
		stencil_end[1] = 2;
		stencil_end[2] = 1;
	}
	
	~OperatorDivergenceLayer() {}
	
	template <typename Lab, typename BlockType>
	void operator()(Lab & lab, const BlockInfo& info, BlockType& o) const
	{
		const Real factor = .5/info.h_gridpoint;
		const int bx = info.index[0]*FluidBlock::sizeX;
		const int by = info.index[1]*FluidBlock::sizeY;
		
		for(int iy=0; iy<FluidBlock::sizeY; ++iy)
			for(int ix=0; ix<FluidBlock::sizeX; ++ix)
			{
				const Real uW = lab(ix-1,iy).u;
				const Real uE = lab(ix+1,iy).u;
				const Real vS = lab(ix,iy-1).v;
				const Real vN = lab(ix,iy+1).v;
				
				divergence(bx + ix, by + iy) = factor * (uE-uW + vN-vS);
			}
	}
};

class OperatorDivergenceSplitLayer : public GenericLabOperator
{
private:
	Layer & divergence;
	const Real rho0, dt;
	const int step;
	
	inline void _mean(const Real c, const Real e, const Real w, const Real n, const Real s, Real& avgE, Real& avgW, Real& avgN, Real& avgS) const
	{
		avgE = .5 * (c + e);
		avgW = .5 * (c + w);
		avgN = .5 * (c + n);
		avgS = .5 * (c + s);
	}
	
public:
	OperatorDivergenceSplitLayer(Layer & divergence, const Real rho0, const Real dt, const int step) : divergence(divergence), rho0(rho0), dt(dt), step(step)
	{
		stencil_start[0] = -1;
		stencil_start[1] = -1;
		stencil_start[2] = 0;
		stencil_end[0] = 2;
		stencil_end[1] = 2;
		stencil_end[2] = 1;
	}
	~OperatorDivergenceSplitLayer() {}
	
	template <typename Lab, typename BlockType>
	void operator()(Lab & lab, const BlockInfo& info, BlockType& o) const
	{
		const Real invH2 = 1./(info.h_gridpoint*info.h_gridpoint);
		const Real factor = rho0 * 0.5/(info.h_gridpoint * dt);
		
		const int bx = info.index[0]*FluidBlock::sizeX;
		const int by = info.index[1]*FluidBlock::sizeY;
		
		for(int iy=0; iy<FluidBlock::sizeY; ++iy)
			for(int ix=0; ix<FluidBlock::sizeX; ++ix)
			{
				FluidElement& phi  = lab(ix  ,iy  );
				FluidElement& phiN = lab(ix  ,iy+1);
				FluidElement& phiS = lab(ix  ,iy-1);
				FluidElement& phiW = lab(ix-1,iy  );
				FluidElement& phiE = lab(ix+1,iy  );
				
				
				Real rhoE, rhoW, rhoN, rhoS;
				_mean(phi.rho, phiE.rho, phiW.rho, phiN.rho, phiS.rho, rhoE, rhoW, rhoN, rhoS);
				
				Real p, pN, pS, pW, pE;
				if (step>=2)
				{
					// p*
					pN = 2*phiN.p - phiN.pOld;
					pS = 2*phiS.p - phiS.pOld;
					pW = 2*phiW.p - phiW.pOld;
					pE = 2*phiE.p - phiE.pOld;
					p  = 2*phi.p  - phi.pOld;
				}
				else
				{
					// pN
					pN = phiN.p;
					pS = phiS.p;
					pW = phiW.p;
					pE = phiE.p;
					p  = phi.p;
				}
				
				Real fN = 1-rho0/rhoN;
				Real fS = 1-rho0/rhoS;
				Real fW = 1-rho0/rhoW;
				Real fE = 1-rho0/rhoE;
				
				divergence(bx + ix, by + iy) = factor * (phiE.u-phiW.u + phiN.v-phiS.v) + invH2 * (fW*pW + fE*pE + fN*pN + fS*pS - (fW+fE+fN+fS)*p);
			}
	}
};

class OperatorDivergence : public GenericLabOperator
{
private:
	double dt;
	
public:
	OperatorDivergence(double dt) : dt(dt)
	{
		stencil_start[0] = -1;
		stencil_start[1] = -1;
		stencil_start[2] = 0;
		stencil_end[0] = 2;
		stencil_end[1] = 2;
		stencil_end[2] = 1;
	}
	~OperatorDivergence() {}
	
	template <typename Lab, typename BlockType>
	void operator()(Lab & lab, const BlockInfo& info, BlockType& o) const
	{
		const Real factor = 0.5/(info.h_gridpoint * dt);
		
		for(int iy=0; iy<FluidBlock::sizeY; ++iy)
			for(int ix=0; ix<FluidBlock::sizeX; ++ix)
			{
				const Real vN = lab(ix  ,iy+1).v;
				const Real vS = lab(ix  ,iy-1).v;
				const Real uW = lab(ix-1,iy  ).u;
				const Real uE = lab(ix+1,iy  ).u;
				o(ix, iy).divU = factor * (uE-uW + vN-vS);
				o(ix, iy).tmp = factor * (uE-uW + vN-vS);
			}
	}
};

class OperatorDivergenceHighOrder : public GenericLabOperator
{
private:
	double dt, rho;

public:
	OperatorDivergenceHighOrder(double dt, double rho) : rho(rho), dt(dt)
	{
		stencil_start[0] = -2;
		stencil_start[1] = -2;
		stencil_start[2] = 0;
		stencil_end[0] = 3;
		stencil_end[1] = 3;
		stencil_end[2] = 1;
	}
	
	~OperatorDivergenceHighOrder() {}
	
	template <typename Lab, typename BlockType>
	void operator()(Lab & lab, const BlockInfo& info, BlockType& o) const
	{
		const Real factor = rho / (12 * info.h_gridpoint * dt);
		
		for(int iy=0; iy<FluidBlock::sizeY; ++iy)
		for(int ix=0; ix<FluidBlock::sizeX; ++ix)
		{
			FluidElement& phi  = lab(ix  ,iy  );
			FluidElement& phiN = lab(ix  ,iy+1);
			FluidElement& phiS = lab(ix  ,iy-1);
			FluidElement& phiW = lab(ix+1,iy  );
			FluidElement& phiE = lab(ix-1,iy  );
			FluidElement& phiN2 = lab(ix  ,iy+2);
			FluidElement& phiS2 = lab(ix  ,iy-2);
			FluidElement& phiW2 = lab(ix+2,iy  );
			FluidElement& phiE2 = lab(ix-2,iy  );
			
			o(ix, iy).divU = factor * (-phiE2.u + 8*phiE.u - 8*phiW.u + phiW2.u - phiN2.v + 8*phiN.v - 8*phiS.v + phiS2.v);
		}
	}
};

class OperatorDivergenceSplit : public GenericLabOperator
{
private:
	double dt, rho0;
	int step;
	
	inline void _mean(const Real c, const Real e, const Real w, const Real n, const Real s, Real& avgE, Real& avgW, Real& avgN, Real& avgS) const
	{
		avgE = .5 * (c + e);
		avgW = .5 * (c + w);
		avgN = .5 * (c + n);
		avgS = .5 * (c + s);
	}
	
	inline void _harmonicAvg(const Real c, const Real e, const Real w, const Real n, const Real s, Real& avgE, Real& avgW, Real& avgN, Real& avgS) const
	{
		avgE = 2. * c * e / (c + e);
		avgW = 2. * c * w / (c + w);
		avgN = 2. * c * n / (c + n);
		avgS = 2. * c * s / (c + s);
	}
	
public:
	OperatorDivergenceSplit(double dt, double rho0, int step) : rho0(rho0), dt(dt), step(step)
	{
		stencil_start[0] = -1;
		stencil_start[1] = -1;
		stencil_start[2] = 0;
		stencil_end[0] = 2;
		stencil_end[1] = 2;
		stencil_end[2] = 1;
	}
	
	~OperatorDivergenceSplit() {}
	
	template <typename Lab, typename BlockType>
	void operator()(Lab & lab, const BlockInfo& info, BlockType& o) const
	{
		const Real invH2 = 1./(info.h_gridpoint*info.h_gridpoint);
		const Real factor = rho0 * 0.5/(info.h_gridpoint * dt);
		
		for(int iy=0; iy<FluidBlock::sizeY; ++iy)
		for(int ix=0; ix<FluidBlock::sizeX; ++ix)
		{
			FluidElement& phi  = lab(ix  ,iy  );
			FluidElement& phiN = lab(ix  ,iy+1);
			FluidElement& phiS = lab(ix  ,iy-1);
			FluidElement& phiW = lab(ix-1,iy  );
			FluidElement& phiE = lab(ix+1,iy  );
			
			
			Real rhoE, rhoW, rhoN, rhoS;
			//_mean(phi.rho, phiE.rho, phiW.rho, phiN.rho, phiS.rho, rhoE, rhoW, rhoN, rhoS);
			_harmonicAvg(phi.rho, phiE.rho, phiW.rho, phiN.rho, phiS.rho, rhoE, rhoW, rhoN, rhoS);
			
			Real p, pN, pS, pW, pE;
			if (step>=2)
			{
				// p*
				pN = 2*phiN.p - phiN.pOld;
				pS = 2*phiS.p - phiS.pOld;
				pW = 2*phiW.p - phiW.pOld;
				pE = 2*phiE.p - phiE.pOld;
				p  = 2*phi.p  - phi.pOld;
			}
			else
			{
				// pN
				pN = phiN.p;
				pS = phiS.p;
				pW = phiW.p;
				pE = phiE.p;
				p  = phi.p;
			}
			
			Real fN = 1-rho0/rhoN;
			Real fS = 1-rho0/rhoS;
			Real fW = 1-rho0/rhoW;
			Real fE = 1-rho0/rhoE;
			
			assert(fN<=1);
			assert(fS<=1);
			assert(fW<=1);
			assert(fE<=1);
			
			o(ix, iy).divU = factor * (phiE.u-phiW.u + phiN.v-phiS.v) + invH2 * (fW*pW + fE*pE + fN*pN + fS*pS - (fW+fE+fN+fS)*p);
			o(ix, iy).tmp = factor * (phiE.u-phiW.u + phiN.v-phiS.v) + invH2 * (fW*pW + fE*pE + fN*pN + fS*pS - (fW+fE+fN+fS)*p);
		}
	}
};

class OperatorDivergenceSplitHighOrder : public GenericLabOperator
{
private:
	double dt, rho0;
	int step;
	
	inline void _mean(const Real c, const Real e, const Real w, const Real n, const Real s, Real& avgE, Real& avgW, Real& avgN, Real& avgS) const
	{
		avgE = .5 * (c + e);
		avgW = .5 * (c + w);
		avgN = .5 * (c + n);
		avgS = .5 * (c + s);
	}
	
	inline void _harmonicAvg(const Real c, const Real e, const Real w, const Real n, const Real s, Real& avgE, Real& avgW, Real& avgN, Real& avgS) const
	{
		avgE = 2. * c * e / (c + e);
		avgW = 2. * c * w / (c + w);
		avgN = 2. * c * n / (c + n);
		avgS = 2. * c * s / (c + s);
	}
	
public:
	OperatorDivergenceSplitHighOrder(double dt, double rho0, int step) : rho0(rho0), dt(dt), step(step)
	{
		stencil_start[0] = -2;
		stencil_start[1] = -2;
		stencil_start[2] = 0;
		stencil_end[0] = 3;
		stencil_end[1] = 3;
		stencil_end[2] = 1;
	}
	
	~OperatorDivergenceSplitHighOrder() {}
	
	template <typename Lab, typename BlockType>
	void operator()(Lab & lab, const BlockInfo& info, BlockType& o) const
	{
		const Real invH2 = 1./(info.h_gridpoint*info.h_gridpoint);
		const Real factor = rho0 / (12 * info.h_gridpoint * dt);
		
		for(int iy=0; iy<FluidBlock::sizeY; ++iy)
		for(int ix=0; ix<FluidBlock::sizeX; ++ix)
		{
			FluidElement& phi  = lab(ix  ,iy  );
			FluidElement& phiN = lab(ix  ,iy+1);
			FluidElement& phiS = lab(ix  ,iy-1);
			FluidElement& phiW = lab(ix+1,iy  );
			FluidElement& phiE = lab(ix-1,iy  );
			FluidElement& phiN2 = lab(ix  ,iy+2);
			FluidElement& phiS2 = lab(ix  ,iy-2);
			FluidElement& phiW2 = lab(ix+2,iy  );
			FluidElement& phiE2 = lab(ix-2,iy  );
			
			o(ix, iy).divU = factor * (-phiE2.u + 8*phiE.u - 8*phiW.u + phiW2.u - phiN2.v + 8*phiN.v - 8*phiS.v + phiS2.v);
			
			Real rhoE, rhoW, rhoN, rhoS;
			//_mean(phi.rho, phiE.rho, phiW.rho, phiN.rho, phiS.rho, rhoE, rhoW, rhoN, rhoS);
			_harmonicAvg(phi.rho, phiE.rho, phiW.rho, phiN.rho, phiS.rho, rhoE, rhoW, rhoN, rhoS);
			
			Real p, pN, pS, pW, pE;
			if (step>=2)
			{
				// p*
				pN = 2*phiN.p - phiN.pOld;
				pS = 2*phiS.p - phiS.pOld;
				pW = 2*phiW.p - phiW.pOld;
				pE = 2*phiE.p - phiE.pOld;
				p  = 2*phi.p  - phi.pOld;
			}
			else
			{
				// pN
				pN = phiN.p;
				pS = phiS.p;
				pW = phiW.p;
				pE = phiE.p;
				p  = phi.p;
			}
			Real fN = 1-rho0/rhoN;
			Real fS = 1-rho0/rhoS;
			Real fW = 1-rho0/rhoW;
			Real fE = 1-rho0/rhoE;
			
			o(ix,iy).divU += invH2 * (fW*pW + fE*pE + fN*pN + fS*pS - (fW+fE+fN+fS)*p);
		}
	}
};

#endif
