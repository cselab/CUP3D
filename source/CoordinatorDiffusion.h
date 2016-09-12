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
	const double mu;
	double dt;
	const int stage;

public:
	OperatorDiffusion(double dt, double mu, const int stage) : mu(mu), dt(dt), stage(stage)
	{
		stencil = StencilInfo(-1,-1,-1, 2,2,2, false, 8, 0,1,2,3,7,8,9,10);

		stencil_start[0] = -1;
		stencil_start[1] = -1;
		stencil_start[2] = -1;

		stencil_end[0] = 2;
		stencil_end[1] = 2;
		stencil_end[2] = 2;

#ifdef _RK2_
		dt = (stage==0) ? dt*.5 : dt;
#endif
	}

	~OperatorDiffusion() {}

	template <typename Lab, typename BlockType>
	void operator()(Lab & lab, const BlockInfo& info, BlockType& o) const
	{
		const double fac = mu * dt / (info.h_gridpoint*info.h_gridpoint);

		// stage 1 of RK2
		if (stage==0)
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

			o(ix,iy,iz).tmpU = phi.u + fac * (phiN.u + phiS.u + phiE.u + phiW.u + phiF.u + phiB.u - phi.u*6.);
			o(ix,iy,iz).tmpV = phi.v + fac * (phiN.v + phiS.v + phiE.v + phiW.v + phiF.v + phiB.v - phi.v*6.);
			o(ix,iy,iz).tmpW = phi.w + fac * (phiN.w + phiS.w + phiE.w + phiW.w + phiF.w + phiB.w - phi.w*6.);
		}

#ifdef _RK2_
		// stage 2 of RK2
		else if (stage==1)
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

			o(ix,iy,iz).tmpU = phi.u + fac * (phiN.tmpU + phiS.tmpU + phiE.tmpU + phiW.tmpU + phiF.tmpU + phiB.tmpU - phi.tmpU*6.);
			o(ix,iy,iz).tmpV = phi.v + fac * (phiN.tmpV + phiS.tmpV + phiE.tmpV + phiW.tmpV + phiF.tmpV + phiB.tmpV - phi.tmpV*6.);
			o(ix,iy,iz).tmpW = phi.w + fac * (phiN.tmpW + phiS.tmpW + phiE.tmpW + phiW.tmpW + phiF.tmpW + phiB.tmpW - phi.tmpW*6.);
		}
#endif // _RK2_
	}
};

template <typename Lab>
class CoordinatorDiffusion : public GenericCoordinator
{
protected:
	const double coeff;
	
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
	
	inline void diffuse(const double dt, const int stage)
	{
		OperatorDiffusion kernel(dt, coeff, stage);
		compute(kernel);
	}
	
public:
	CoordinatorDiffusion(const double coeff, FluidGridMPI * grid) : GenericCoordinator(grid), coeff(coeff)
	{
	}
	
	void operator()(const double dt)
	{
		check("diffusion - start");
		
		diffuse(dt,0);
#ifdef _RK2_
		diffuse(dt,1);
#endif
		update();
		
		check("diffusion - end");
	}
	
	string getName()
	{
		return "Diffusion";
	}
};


class CoordinatorDiffusionTimeTest : public GenericCoordinator
{
protected:
	const double coeff;
	double time;
	const double freq;
	
	double _analytical(double px, double py, double pz, double t)
	{
		return sin(px*2.*freq*M_PI) * sin(py*2.*freq*M_PI) * sin(pz*2.*freq*M_PI) * exp(-4.*3*freq*freq*coeff*M_PI*M_PI*t);
	}
	
	double _analyticalRHS(double px, double py, double pz, double t)
	{
		return -freq*freq*4.*3.*M_PI*M_PI*_analytical(px,py,pz,t);
	}
	
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
	
	inline void update()
	{
		const int N = vInfo.size();
		
#pragma omp parallel for schedule(static)
		for(int i=0; i<N; i++)
		{
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
	
	inline void diffuse(const double dt, const int stage)
	{
		const int N = vInfo.size();
		
#pragma omp parallel for schedule(static)
		for(int i=0; i<N; i++)
		{
			BlockInfo info = vInfo[i];
			FluidBlock& b = *(FluidBlock*)info.ptrBlock;
			const double prefactor = coeff * dt * ((stage==0)?.5:1);
			
			for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
				for(int iy=0; iy<FluidBlock::sizeY; ++iy)
					for(int ix=0; ix<FluidBlock::sizeX; ++ix)
					{
						double p[3];
						info.pos(p, ix, iy, iz);
						
						if (stage==0)
							b(ix,iy,iz).tmpU = b(ix,iy,iz).u + prefactor * _analyticalRHS(p[0],p[1],p[2],time);
						else if (stage==1)
							b(ix,iy,iz).tmpU = b(ix,iy,iz).u + prefactor * _analyticalRHS(p[0],p[1],p[2],time+dt*.5);
						b(ix,iy,iz).tmpV = b(ix,iy,iz).v;
						b(ix,iy,iz).tmpW = b(ix,iy,iz).w;
					}
		}
	}
	
public:
	CoordinatorDiffusionTimeTest(const double coeff, const double freq, FluidGridMPI * grid) : GenericCoordinator(grid), coeff(coeff), freq(freq), time(0)
	{
	}
	
	void operator()(const double dt)
	{
		check("diffusion - start");
		
		reset();
		diffuse(dt,0);
#ifdef _RK2_
		diffuse(dt,1);
#endif
		update();
		time+=dt;
		
		check("diffusion - end");
	}
	
	string getName()
	{
		return "DiffusionTimeTest";
	}
};
#endif
