//
//  OperatorGradP.h
//  CubismUP_3D
//
//	Operates on
//		u, v, tmp (diagnostic)
//
//  Created by Christian Conti on 1/15/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef CubismUP_3D_OperatorGradP_h
#define CubismUP_3D_OperatorGradP_h

#include "GenericOperator.h"

class OperatorGradP : public GenericLabOperator
{
private:
	double dt;
	
public:
	OperatorGradP(double dt) : dt(dt)
	{
		stencil_start[0] = -1;
		stencil_start[1] = -1;
		stencil_start[2] = 0;
		stencil_end[0] = 2;
		stencil_end[1] = 2;
		stencil_end[2] = 1;
	}
	
	~OperatorGradP() {}
	
	template <typename Lab, typename BlockType>
	void operator()(Lab & lab, const BlockInfo& info, BlockType& o) const
	{
		const double prefactor = -.5 * dt / (info.h_gridpoint);
		
		for(int iy=0; iy<FluidBlock::sizeY; ++iy)
			for(int ix=0; ix<FluidBlock::sizeX; ++ix)
			{
				const Real divUN = lab(ix  ,iy+1).divU;
				const Real divUS = lab(ix  ,iy-1).divU;
				const Real divUE = lab(ix+1,iy  ).divU;
				const Real divUW = lab(ix-1,iy  ).divU;
				
				// divU contains the pressure correction after the Poisson solver
				o(ix,iy).u += prefactor * (divUE - divUW) / lab(ix,iy).rho;
				o(ix,iy).v += prefactor * (divUN - divUS) / lab(ix,iy).rho;
				
				assert(lab(ix,iy).rho > 0);
				assert(!std::isnan(o(ix,iy).u));
				assert(!std::isnan(o(ix,iy).v));
			}
	}
};

class OperatorGradPSplit : public GenericLabOperator
{
private:
	double rho0;
	double dt;
	int step;
	
public:
	OperatorGradPSplit(double dt, double rho0, int step) : rho0(rho0), dt(dt), step(step)
	{
		stencil_start[0] = -1;
		stencil_start[1] = -1;
		stencil_start[2] = 0;
		stencil_end[0] = 2;
		stencil_end[1] = 2;
		stencil_end[2] = 1;
	}
	
	~OperatorGradPSplit() {}
	
	template <typename Lab, typename BlockType>
	void operator()(Lab & lab, const BlockInfo& info, BlockType& o) const
	{
		const double dh = info.h_gridpoint;
		const double prefactor = -.5 * dt / dh;
		
		for(int iy=0; iy<FluidBlock::sizeY; ++iy)
			for(int ix=0; ix<FluidBlock::sizeX; ++ix)
			{
				const FluidElement& phi  = lab(ix  ,iy  );
				const FluidElement& phiN = lab(ix  ,iy+1);
				const FluidElement& phiS = lab(ix  ,iy-1);
				const FluidElement& phiE = lab(ix+1,iy  );
				const FluidElement& phiW = lab(ix-1,iy  );
				
				Real pN, pS, pW, pE;
				
				if (step>=2)
				{
					// p*
					pN = 2*phiN.p - phiN.pOld;
					pS = 2*phiS.p - phiS.pOld;
					pW = 2*phiW.p - phiW.pOld;
					pE = 2*phiE.p - phiE.pOld;
				}
				else
				{
					// pN
					pN = phiN.p;
					pS = phiS.p;
					pW = phiW.p;
					pE = phiE.p;
				}
				
				// divU contains the pressure correction after the Poisson solver
				o(ix,iy).u += prefactor/rho0 * (phiE.divU - phiW.divU);
				o(ix,iy).v += prefactor/rho0 * (phiN.divU - phiS.divU);
				
				// add the split explicit term
				o(ix,iy).u += prefactor * (pE - pW) * (1./phi.rho - 1/rho0);
				o(ix,iy).v += prefactor * (pN - pS) * (1./phi.rho - 1/rho0);
				
				o(ix,iy).tmp = prefactor/rho0 * (phiN.divU - phiS.divU) + prefactor * (pN - pS) * (1./phi.rho - 1/rho0);
			}
	}
};

class OperatorGradPSplitHighOrder : public GenericLabOperator
{
private:
	double rho0;
	double dt;
	int step;
	
public:
	OperatorGradPSplitHighOrder(double dt, double rho0, int step) : rho0(rho0), dt(dt), step(step)
	{
		stencil_start[0] = -2;
		stencil_start[1] = -2;
		stencil_start[2] = 0;
		stencil_end[0] = 3;
		stencil_end[1] = 3;
		stencil_end[2] = 1;
	}
	
	~OperatorGradPSplitHighOrder() {}
	
	template <typename Lab, typename BlockType>
	void operator()(Lab & lab, const BlockInfo& info, BlockType& o) const
	{
		const double dh = info.h_gridpoint;
		const double prefactor = -dt / (12.*dh);
		
		for(int iy=0; iy<FluidBlock::sizeY; ++iy)
			for(int ix=0; ix<FluidBlock::sizeX; ++ix)
			{
				const FluidElement& phi  = lab(ix  ,iy  );
				const FluidElement& phiN = lab(ix  ,iy+1);
				const FluidElement& phiS = lab(ix  ,iy-1);
				const FluidElement& phiE = lab(ix+1,iy  );
				const FluidElement& phiW = lab(ix-1,iy  );
				const FluidElement& phiN2 = lab(ix  ,iy+2);
				const FluidElement& phiS2 = lab(ix  ,iy-2);
				const FluidElement& phiE2 = lab(ix+2,iy  );
				const FluidElement& phiW2 = lab(ix-2,iy  );
				
				Real p, pN, pS, pW, pE, pN2, pS2, pW2, pE2;
				if (step>=2)
				{
					// p*
					pN = 2*phiN.p - phiN.pOld;
					pS = 2*phiS.p - phiS.pOld;
					pW = 2*phiW.p - phiW.pOld;
					pE = 2*phiE.p - phiE.pOld;
					pN2 = 2*phiN2.p - phiN2.pOld;
					pS2 = 2*phiS2.p - phiS2.pOld;
					pW2 = 2*phiW2.p - phiW2.pOld;
					pE2 = 2*phiE2.p - phiE2.pOld;
				}
				else
				{
					// pN
					pN = phiN.p;
					pS = phiS.p;
					pW = phiW.p;
					pE = phiE.p;
					pN2 = phiN2.p;
					pS2 = phiS2.p;
					pW2 = phiW2.p;
					pE2 = phiE2.p;
				}
				
				// divU contains the pressure correction after the Poisson solver
				o(ix,iy).u += prefactor/rho0 * (-phiE2.divU + 8*phiE.divU - 8*phiW.divU + phiW2.divU);
				o(ix,iy).v += prefactor/rho0 * (-phiN2.divU + 8*phiN.divU - 8*phiS.divU + phiS2.divU);
				
				// add the split explicit term
				o(ix,iy).u += prefactor * (-pE2 + 8*pE - 8*pW + pW2) * (1./phi.rho - 1/rho0);
				o(ix,iy).v += prefactor * (-pN2 + 8*pN - 8*pS + pS2) * (1./phi.rho - 1/rho0);
			}
	}
};

#endif
