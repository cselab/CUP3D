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
		stencil_start[2] = -1;
		stencil_end[0] = 2;
		stencil_end[1] = 2;
		stencil_end[2] = 2;
	}
	
	~OperatorDivergenceLayer() {}
	
	template <typename Lab, typename BlockType>
	void operator()(Lab & lab, const BlockInfo& info, BlockType& o) const
	{
		const Real factor = .5/info.h_gridpoint;
		const int bx = info.index[0]*FluidBlock::sizeX;
		const int by = info.index[1]*FluidBlock::sizeY;
		const int bz = info.index[2]*FluidBlock::sizeZ;
		
		for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
			for(int iy=0; iy<FluidBlock::sizeY; ++iy)
				for(int ix=0; ix<FluidBlock::sizeX; ++ix)
				{
					const Real uW = lab(ix-1,iy,iz).u;
					const Real uE = lab(ix+1,iy,iz).u;
					const Real vS = lab(ix,iy-1,iz).v;
					const Real vN = lab(ix,iy+1,iz).v;
					const Real wF = lab(ix,iy,iz-1).w;
					const Real wB = lab(ix,iy,iz+1).w;
					
					divergence(bx + ix, by + iy) = factor * (uE-uW + vN-vS + wB-wF);
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
		stencil_start[2] = -1;
		stencil_end[0] = 2;
		stencil_end[1] = 2;
		stencil_end[2] = 2;
	}
	~OperatorDivergence() {}
	
	template <typename Lab, typename BlockType>
	void operator()(Lab & lab, const BlockInfo& info, BlockType& o) const
	{
		const Real factor = 0.5/(info.h_gridpoint * dt);
		
		for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
			for(int iy=0; iy<FluidBlock::sizeY; ++iy)
				for(int ix=0; ix<FluidBlock::sizeX; ++ix)
				{
					const Real uW = lab(ix-1,iy  ,iz  ).u;
					const Real uE = lab(ix+1,iy  ,iz  ).u;
					const Real vS = lab(ix  ,iy-1,iz  ).v;
					const Real vN = lab(ix  ,iy+1,iz  ).v;
					const Real wF = lab(ix  ,iy  ,iz-1).w;
					const Real wB = lab(ix  ,iy  ,iz+1).w;
					o(ix, iy, iz).divU = factor * (uE-uW + vN-vS + wB-wF);
					o(ix, iy, iz).tmp  = factor * (uE-uW + vN-vS + wB-wF);
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
		stencil_start[2] = -2;
		stencil_end[0] = 3;
		stencil_end[1] = 3;
		stencil_end[2] = 3;
	}
	
	~OperatorDivergenceHighOrder() {}
	
	template <typename Lab, typename BlockType>
	void operator()(Lab & lab, const BlockInfo& info, BlockType& o) const
	{
		const Real factor = rho / (12 * info.h_gridpoint * dt);
		
		for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
			for(int iy=0; iy<FluidBlock::sizeY; ++iy)
				for(int ix=0; ix<FluidBlock::sizeX; ++ix)
				{
					FluidElement& phi  = lab(ix  ,iy  ,iz  );
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
					
					o(ix, iy, iz).divU = factor * (-phiE2.u + 8*phiE.u - 8*phiW.u + phiW2.u - phiN2.v + 8*phiN.v - 8*phiS.v + phiS2.v - phiB2.w + 8*phiB.w - 8*phiF.w + phiF2.w);
				}
	}
};

class OperatorDivergenceSplit : public GenericLabOperator
{
private:
	double dt, rho0;
	int step;
	
	inline void _mean(const Real c, const Real w, const Real e, const Real s, const Real n, const Real f, const Real b, Real& avgW, Real& avgE, Real& avgS, Real& avgN, Real& avgF, Real& avgB) const
	{
		avgW = .5 * (c + w);
		avgE = .5 * (c + e);
		avgS = .5 * (c + s);
		avgN = .5 * (c + n);
		avgF = .5 * (c + f);
		avgB = .5 * (c + b);
	}
	
	inline void _harmonicAvg(const Real c, const Real w, const Real e, const Real s, const Real n, const Real f, const Real b, Real& avgW, Real& avgE, Real& avgS, Real& avgN, Real& avgF, Real& avgB) const
	{
		avgW = 2. * c * w / (c + w);
		avgE = 2. * c * e / (c + e);
		avgS = 2. * c * s / (c + s);
		avgN = 2. * c * n / (c + n);
		avgF = 2. * c * f / (c + f);
		avgB = 2. * c * b / (c + b);
	}
	
public:
	OperatorDivergenceSplit(double dt, double rho0, int step) : rho0(rho0), dt(dt), step(step)
	{
		stencil = StencilInfo(-1,-1,-1, 2,2,2, false, 6, 0,1,2,3,5,6);
		
		stencil_start[0] = -1;
		stencil_start[1] = -1;
		stencil_start[2] = -1;
		stencil_end[0] = 2;
		stencil_end[1] = 2;
		stencil_end[2] = 2;
	}
	
	~OperatorDivergenceSplit() {}
	
	template <typename Lab, typename BlockType>
	void operator()(Lab & lab, const BlockInfo& info, BlockType& o) const
	{
		const Real invH2 = 1./(info.h_gridpoint*info.h_gridpoint);
		const Real factor = rho0 * 0.5/(info.h_gridpoint * dt);
		
		for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
			for(int iy=0; iy<FluidBlock::sizeY; ++iy)
				for(int ix=0; ix<FluidBlock::sizeX; ++ix)
				{
					FluidElement& phi  = lab(ix  ,iy  ,iz  );
					FluidElement& phiW = lab(ix-1,iy  ,iz  );
					FluidElement& phiE = lab(ix+1,iy  ,iz  );
					FluidElement& phiS = lab(ix  ,iy-1,iz  );
					FluidElement& phiN = lab(ix  ,iy+1,iz  );
					FluidElement& phiF = lab(ix  ,iy  ,iz-1);
					FluidElement& phiB = lab(ix  ,iy  ,iz+1);
					
					
					Real rhoW, rhoE, rhoS, rhoN, rhoF, rhoB;
					//_mean(phi.rho, phiE.rho, phiW.rho, phiN.rho, phiS.rho, rhoE, rhoW, rhoN, rhoS);
					_harmonicAvg(phi.rho, phiW.rho, phiE.rho, phiS.rho, phiN.rho, phiF.rho, phiB.rho, rhoW, rhoE, rhoS, rhoN, rhoF, rhoB);
					
					Real p, pW, pE, pS, pN, pF, pB;
					if (step>=2)
					{
						// p*
						pW = 2*phiW.p - phiW.pOld;
						pE = 2*phiE.p - phiE.pOld;
						pS = 2*phiS.p - phiS.pOld;
						pN = 2*phiN.p - phiN.pOld;
						pF = 2*phiF.p - phiF.pOld;
						pB = 2*phiB.p - phiB.pOld;
						p  = 2*phi.p  - phi.pOld;
					}
					else
					{
						// pN
						pW = phiW.p;
						pE = phiE.p;
						pS = phiS.p;
						pN = phiN.p;
						pF = phiF.p;
						pB = phiB.p;
						p  = phi.p;
					}
					
					Real fW = 1-rho0/rhoW;
					Real fE = 1-rho0/rhoE;
					Real fS = 1-rho0/rhoS;
					Real fN = 1-rho0/rhoN;
					Real fF = 1-rho0/rhoF;
					Real fB = 1-rho0/rhoB;
					
					assert(fW<=1);
					assert(fE<=1);
					assert(fS<=1);
					assert(fN<=1);
					assert(fF<=1);
					assert(fB<=1);
					
					o(ix, iy, iz).divU = factor * (phiE.u-phiW.u + phiN.v-phiS.v + phiB.w-phiF.w) + invH2 * (fW*pW + fE*pE + fS*pS + fN*pN + fF*pF + fB*pB - (fW+fE+fS+fN+fF+fB)*p);
					o(ix, iy, iz).tmp  = factor * (phiE.u-phiW.u + phiN.v-phiS.v + phiB.w-phiF.w) + invH2 * (fW*pW + fE*pE + fS*pS + fN*pN + fF*pF + fB*pB - (fW+fE+fS+fN+fF+fB)*p);
				}
	}
};

class OperatorDivergenceSplitFFTW : public GenericLabOperator
{
private:
	double dt, rho0;
	int step;
	Real * data; // FFTW data structure
	int dim[3];
	
	inline void _mean(const Real c, const Real w, const Real e, const Real s, const Real n, const Real f, const Real b, Real& avgW, Real& avgE, Real& avgS, Real& avgN, Real& avgF, Real& avgB) const
	{
		avgW = .5 * (c + w);
		avgE = .5 * (c + e);
		avgS = .5 * (c + s);
		avgN = .5 * (c + n);
		avgF = .5 * (c + f);
		avgB = .5 * (c + b);
	}
	
	inline void _harmonicAvg(const Real c, const Real w, const Real e, const Real s, const Real n, const Real f, const Real b, Real& avgW, Real& avgE, Real& avgS, Real& avgN, Real& avgF, Real& avgB) const
	{
		avgW = 2. * c * w / (c + w);
		avgE = 2. * c * e / (c + e);
		avgS = 2. * c * s / (c + s);
		avgN = 2. * c * n / (c + n);
		avgF = 2. * c * f / (c + f);
		avgB = 2. * c * b / (c + b);
	}
	
public:
	
	OperatorDivergenceSplitFFTW(double dt, double rho0, int step, Real * data, int dim[3]) : rho0(rho0), dt(dt), step(step), data(data), dim{dim[0],dim[1],dim[2]}
	{
		stencil = StencilInfo(-1,-1,-1, 2,2,2, false, 6, 0,1,2,3,5,6);
		
		stencil_start[0] = -1;
		stencil_start[1] = -1;
		stencil_start[2] = -1;
		stencil_end[0] = 2;
		stencil_end[1] = 2;
		stencil_end[2] = 2;
	}
	
	~OperatorDivergenceSplitFFTW() {}
	
	template <typename Lab, typename BlockType>
	void operator()(Lab & lab, const BlockInfo& info, BlockType& o)
	{
		const Real invH2 = 1./(info.h_gridpoint*info.h_gridpoint);
		const Real factor = rho0 * 0.5/(info.h_gridpoint * dt);
		
		const size_t offset = FluidBlock::sizeZ*info.index[2] + dim[2]*FluidBlock::sizeY*info.index[1] + dim[2]*dim[1]*FluidBlock::sizeX*info.index[0];
		
		for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
			for(int iy=0; iy<FluidBlock::sizeY; ++iy)
				for(int ix=0; ix<FluidBlock::sizeX; ++ix)
				{
					FluidElement& phi  = lab(ix  ,iy  ,iz  );
					FluidElement& phiW = lab(ix-1,iy  ,iz  );
					FluidElement& phiE = lab(ix+1,iy  ,iz  );
					FluidElement& phiS = lab(ix  ,iy-1,iz  );
					FluidElement& phiN = lab(ix  ,iy+1,iz  );
					FluidElement& phiF = lab(ix  ,iy  ,iz-1);
					FluidElement& phiB = lab(ix  ,iy  ,iz+1);
					
					
					Real rhoW, rhoE, rhoS, rhoN, rhoF, rhoB;
					//_mean(phi.rho, phiE.rho, phiW.rho, phiN.rho, phiS.rho, rhoE, rhoW, rhoN, rhoS);
					_harmonicAvg(phi.rho, phiW.rho, phiE.rho, phiS.rho, phiN.rho, phiF.rho, phiB.rho, rhoW, rhoE, rhoS, rhoN, rhoF, rhoB);
					
					Real p, pW, pE, pS, pN, pF, pB;
					if (step>=2)
					{
						// p*
						pW = 2*phiW.p - phiW.pOld;
						pE = 2*phiE.p - phiE.pOld;
						pS = 2*phiS.p - phiS.pOld;
						pN = 2*phiN.p - phiN.pOld;
						pF = 2*phiF.p - phiF.pOld;
						pB = 2*phiB.p - phiB.pOld;
						p  = 2*phi.p  - phi.pOld;
					}
					else
					{
						// pN
						pW = phiW.p;
						pE = phiE.p;
						pS = phiS.p;
						pN = phiN.p;
						pF = phiF.p;
						pB = phiB.p;
						p  = phi.p;
					}
					
					Real fW = 1-rho0/rhoW;
					Real fE = 1-rho0/rhoE;
					Real fS = 1-rho0/rhoS;
					Real fN = 1-rho0/rhoN;
					Real fF = 1-rho0/rhoF;
					Real fB = 1-rho0/rhoB;
					
					assert(fW<=1);
					assert(fE<=1);
					assert(fS<=1);
					assert(fN<=1);
					assert(fF<=1);
					assert(fB<=1);
					
					const size_t dest_index = offset + iz + dim[2]*iy + dim[1]*dim[2]*ix;
					data[dest_index] = factor * (phiE.u-phiW.u + phiN.v-phiS.v + phiB.w-phiF.w) + invH2 * (fW*pW + fE*pE + fS*pS + fN*pN + fF*pF + fB*pB - (fW+fE+fS+fN+fF+fB)*p);
					
					//o(ix, iy, iz).divU = factor * (phiE.u-phiW.u + phiN.v-phiS.v + phiB.w-phiF.w) + invH2 * (fW*pW + fE*pE + fS*pS + fN*pN + fF*pF + fB*pB - (fW+fE+fS+fN+fF+fB)*p);
					//o(ix, iy, iz).tmp  = factor * (phiE.u-phiW.u + phiN.v-phiS.v + phiB.w-phiF.w) + invH2 * (fW*pW + fE*pE + fS*pS + fN*pN + fF*pF + fB*pB - (fW+fE+fS+fN+fF+fB)*p);
				}
	}
};

class OperatorDivergenceSplitHighOrder : public GenericLabOperator
{
private:
	double dt, rho0;
	int step;
	
	inline void _mean(const Real c, const Real w, const Real e, const Real s, const Real n, const Real f, const Real b, Real& avgW, Real& avgE, Real& avgS, Real& avgN, Real& avgF, Real& avgB) const
	{
		avgW = .5 * (c + w);
		avgE = .5 * (c + e);
		avgS = .5 * (c + s);
		avgN = .5 * (c + n);
		avgF = .5 * (c + f);
		avgB = .5 * (c + b);
	}
	
	inline void _harmonicAvg(const Real c, const Real w, const Real e, const Real s, const Real n, const Real f, const Real b, Real& avgW, Real& avgE, Real& avgS, Real& avgN, Real& avgF, Real& avgB) const
	{
		avgW = 2. * c * w / (c + w);
		avgE = 2. * c * e / (c + e);
		avgS = 2. * c * s / (c + s);
		avgN = 2. * c * n / (c + n);
		avgF = 2. * c * f / (c + f);
		avgB = 2. * c * b / (c + b);
	}
	
public:
	OperatorDivergenceSplitHighOrder(double dt, double rho0, int step) : rho0(rho0), dt(dt), step(step)
	{
		stencil = StencilInfo(-2,-2,-2, 3,3,3, false, 6, 0,1,2,3,5,6);
		
		stencil_start[0] = -2;
		stencil_start[1] = -2;
		stencil_start[2] = -2;
		stencil_end[0] = 3;
		stencil_end[1] = 3;
		stencil_end[2] = 3;
	}
	
	~OperatorDivergenceSplitHighOrder() {}
	
	template <typename Lab, typename BlockType>
	void operator()(Lab & lab, const BlockInfo& info, BlockType& o) const
	{
		const Real invH2 = 1./(info.h_gridpoint*info.h_gridpoint);
		const Real factor = rho0 / (12 * info.h_gridpoint * dt);
		
		for (int iz=0; iz<FluidBlock::sizeY; ++iz)
			for (int iy=0; iy<FluidBlock::sizeY; ++iy)
				for (int ix=0; ix<FluidBlock::sizeX; ++ix)
				{
					FluidElement& phi  = lab(ix  ,iy  ,iz  );
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
					
					o(ix, iy).divU = factor * (-phiE2.u + 8*phiE.u - 8*phiW.u + phiW2.u - phiN2.v + 8*phiN.v - 8*phiS.v + phiS2.v - phiB2.v + 8*phiB.v - 8*phiF.v + phiF2.v);
					
					Real rhoW, rhoE, rhoS, rhoN, rhoF, rhoB;
					//_mean(phi.rho, phiE.rho, phiW.rho, phiN.rho, phiS.rho, rhoE, rhoW, rhoN, rhoS);
					_harmonicAvg(phi.rho, phiW.rho, phiE.rho, phiS.rho, phiN.rho, phiF.rho, phiB.rho, rhoW, rhoE, rhoS, rhoN, rhoF, rhoB);
					
					Real p, pN, pS, pW, pE, pB, pF;
					if (step>=2)
					{
						// p*
						pN = 2*phiN.p - phiN.pOld;
						pS = 2*phiS.p - phiS.pOld;
						pW = 2*phiW.p - phiW.pOld;
						pE = 2*phiE.p - phiE.pOld;
						pB = 2*phiB.p - phiB.pOld;
						pF = 2*phiF.p - phiF.pOld;
						p  = 2*phi.p  - phi.pOld;
					}
					else
					{
						// pN
						pN = phiN.p;
						pS = phiS.p;
						pW = phiW.p;
						pE = phiE.p;
						pB = phiB.p;
						pF = phiF.p;
						p  = phi.p;
					}
					Real fN = 1-rho0/rhoN;
					Real fS = 1-rho0/rhoS;
					Real fW = 1-rho0/rhoW;
					Real fE = 1-rho0/rhoE;
					Real fF = 1-rho0/rhoF;
					Real fB = 1-rho0/rhoB;
					
					o(ix,iy).divU += invH2 * (fW*pW + fE*pE + fN*pN + fS*pS + fF*pF + fB*pB - (fW+fE+fN+fS+fF+fB)*p);
				}
	}
};

#endif
