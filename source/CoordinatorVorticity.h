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
		stencil = StencilInfo(-1,-1,-1, 2,2,2, false, 3, 0,1,2);

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
		const Real inv2h = .5 / info.h_gridpoint;

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

class OperatorDiagnostis : public GenericOperator
{
public:
	array<Real,12>* const quantities;

	OperatorDiagnostis(array<Real,12>* const local) : quantities(local) {}

	~OperatorDiagnostis() {}

	void operator()(const BlockInfo& info, FluidBlock& block) const
	{
		const Real dv = std::pow(info.h_gridpoint,3);
		for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
		for(int iy=0; iy<FluidBlock::sizeY; ++iy)
		for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
			Real x[3];
			info.pos(x, ix, iy, iz);
			const Real w[3] = {
					b(ix,iy,iz).omega[0],
					b(ix,iy,iz).omega[1],
					b(ix,iy,iz).omega[2]
			};
			const Real u[3] = {
					b(ix,iy,iz).u[0],
					b(ix,iy,iz).u[1],
					b(ix,iy,iz).u[2]
			};
			const Real omegasq = w[0]*w[0] + w[1]*w[1] + w[2]*w[2];
			const Real velsq   = u[0]*u[0] + u[1]*u[1] + u[2]*u[2];
			//circulation :
			(*quantities)[0] += w[0];
			(*quantities)[1] += w[1];
			(*quantities)[2] += w[2];
			//linear impulse :
			(*quantities)[3] += x[1]*w[2]-x[2]*w[1];
			(*quantities)[4] += x[2]*w[0]-x[0]*w[2];
			(*quantities)[5] += x[0]*w[1]-x[1]*w[0];
			//angular impulse :
			(*quantities)[6] += x[0]*(x[1]*w[1]+x[2]*w[2]) - (x[1]*x[1]+x[2]*x[2])*w[0];
			(*quantities)[7] += x[1]*(x[0]*w[0]+x[2]*w[2]) - (x[0]*x[0]+x[2]*x[2])*w[1];
			(*quantities)[8] += x[2]*(x[0]*w[0]+x[1]*w[1]) - (x[0]*x[0]+x[1]*x[1])*w[2];
			//helicity
			(*quantities)[9] += (w[0]*u[0] + w[1]*u[1] + w[2]*u[2]);
			//enstrophy
			(*quantities)[10]+= omegasq;
			//maxvor
			(*quantities)[11] = std::max((*quantities)[11],omegasq);
		}

		for(int i=0;i<3;i++) {
			(*quantities)[  i] *= dv;
			(*quantities)[3+i] *= 0.5*dv;
			(*quantities)[6+i] *= dv/3.;
		}
		(*quantities)[9] *= dv;
		(*quantities)[10]*= dv;
	}
	/*
	FILE * f = fopen("diagnostics.dat", "a");
	if(step_id == 0 && !bRESTART)
		fprintf(f,"%s  %s  %s  %s  %s  %s  %s  %s  %s  %s\n",
				"# step_id","time","dt","circ","linImpX","linImpY","angImp","Eng","Ens","Maxvor");

	fprintf(f, "%d  %10.10e  %10.10e  %10.10e  %10.10e  %10.10e  %10.10e  %10.10e  %10.10e  %10.10e\n",
			step_id, t, dt, diags.data.circ,
			diags.data.linImpulse[0], diags.data.linImpulse[1],diags.data.angImpulse,
			keSolver,diags.data.ens,diags.data.maxvor);

	fclose(f);
	*/
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
	
	void operator()(const Real dt)
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
