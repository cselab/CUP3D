//
//  OperatorSplitP.h
//  CubismUP_3D
//
//	Operates on
//		divU
//
//  Created by Christian Conti on 1/28/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef CubismUP_3D_OperatorSplitP_h
#define CubismUP_3D_OperatorSplitP_h

#include "GenericOperator.h"

class OperatorSplitP : public GenericLabOperator
{
private:
	double rho0;
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
	OperatorSplitP(double rho0, int step) : rho0(rho0), step(step)
	{
		stencil_start[0] = -1;
		stencil_start[1] = -1;
		stencil_start[2] = 0;
		stencil_end[0] = 2;
		stencil_end[1] = 2;
		stencil_end[2] = 1;
	}
	
	~OperatorSplitP() {}
	
	template <typename Lab, typename BlockType>
	void operator()(Lab & lab, const BlockInfo& info, BlockType& o) const
	{
		const Real dh = info.h_gridpoint;
		const Real invH2 = 1./(dh*dh);
		
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
				
				Real fN = 1.-rho0/rhoN;
				Real fS = 1.-rho0/rhoS;
				Real fW = 1.-rho0/rhoW;
				Real fE = 1.-rho0/rhoE;
				
				o(ix,iy).divU += invH2 * (fE*pE + fW*pW + fN*pN + fS*pS - (fW+fE+fN+fS)*p);
			}
	}
};

#endif
