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
	void operator()(const Lab & lab, const BlockInfo& info, BlockType& o) const
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

class OperatorDiagnostics : public GenericOperator
{
public:
	vector<array<Real,12>>* const quantities;

	OperatorDiagnostics(vector<array<Real,12>>* const local) : quantities(local) {}

	void operator()(const BlockInfo& info, FluidBlock& b) const
	{
		Real circ[3]   = {0.,0.,0.};
		Real linimp[3] = {0.,0.,0.};
		Real angimp[3] = {0.,0.,0.};
		Real helicity(0), enstrophy(0), maxvortSq(0);

		for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
		for(int iy=0; iy<FluidBlock::sizeY; ++iy)
		for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
			Real x[3];
			info.pos(x, ix, iy, iz);
			const Real w[3] = {
					b(ix,iy,iz).tmpU,
					b(ix,iy,iz).tmpV,
					b(ix,iy,iz).tmpW
			};
			const Real u[3] = {
					b(ix,iy,iz).u,
					b(ix,iy,iz).v,
					b(ix,iy,iz).w
			};
			const Real omegasq = w[0]*w[0] + w[1]*w[1] + w[2]*w[2];
			const Real velsq   = u[0]*u[0] + u[1]*u[1] + u[2]*u[2];
			//circulation :
			circ[0] += w[0];
			circ[1] += w[1];
			circ[2] += w[2];
			//linear impulse :
			linimp[0] += x[1]*w[2]-x[2]*w[1];
			linimp[1] += x[2]*w[0]-x[0]*w[2];
			linimp[2] += x[0]*w[1]-x[1]*w[0];
			//angular impulse :
			angimp[0] += x[0]*(x[1]*w[1]+x[2]*w[2]) - (x[1]*x[1]+x[2]*x[2])*w[0];
			angimp[1] += x[1]*(x[0]*w[0]+x[2]*w[2]) - (x[0]*x[0]+x[2]*x[2])*w[1];
			angimp[2] += x[2]*(x[0]*w[0]+x[1]*w[1]) - (x[0]*x[0]+x[1]*x[1])*w[2];

			helicity += (w[0]*u[0] + w[1]*u[1] + w[2]*u[2]);
			enstrophy+= omegasq;
			maxvortSq = std::max(maxvortSq,omegasq);
		}

		const Real dv = std::pow(info.h_gridpoint,3);
		const int tid = omp_get_thread_num();
		for(int i=0;i<3;i++) {
			(*quantities)[tid][  i] += circ[i]*dv;
			(*quantities)[tid][3+i] += linimp[i]*0.5*dv;
			(*quantities)[tid][6+i] += angimp[i]*dv/3.;
		}
		(*quantities)[tid][9] += helicity*dv;
		(*quantities)[tid][10]+= enstrophy*dv;
		(*quantities)[tid][11] = std::max(maxvortSq,(*quantities)[tid][11]);
	}
};

template <typename Lab>
class CoordinatorVorticity : public GenericCoordinator
{
public:
	CoordinatorVorticity(FluidGridMPI * grid) : GenericCoordinator(grid) { }
	
	void operator()(const Real dt)
	{
		check("vorticity - start");
		OperatorVorticity kernel;
		compute(kernel);
		check("vorticity - end");
	}
	
	string getName() { return "Vorticity"; }
};

class CoordinatorDiagnostics : public GenericCoordinator
{
private:
	const Real t;
	const int step;
	int rank;
public:
	CoordinatorDiagnostics(FluidGridMPI * grid, const Real t, const int step)
	: GenericCoordinator(grid), t(t), step(step)
	{
		MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	}

	void operator()(const Real dt)
	{
	   const int nthreads = omp_get_max_threads();
	   vector<array<Real,12>> partialSums(nthreads);
      const int N = vInfo.size();
#pragma omp parallel
      {
         OperatorDiagnostics kernel(&partialSums);
#pragma omp for schedule(static)
         for(int i=0; i<N; i++) {
            BlockInfo info = vInfo[i];
            FluidBlock& b = *(FluidBlock*)info.ptrBlock;
            kernel(info, b);
         }
      }

		double localSum[11], globalSum[11], maxVortHere(0), maxVort(0);
		for(int i=0; i<nthreads; i++) {
			for(int j=0; j<11; j++)
				localSum[j] += (double)partialSums[i][j];
			maxVortHere = std::max((double)partialSums[i][11],maxVortHere);
		}
		MPI::COMM_WORLD.Allreduce( localSum, globalSum, 11, MPI::DOUBLE, MPI::SUM);
		MPI::COMM_WORLD.Allreduce(&maxVortHere,&maxVort, 1, MPI::DOUBLE, MPI::MAX);

		if(rank==0) {
			FILE * f = fopen("diagnostics.dat", "a");
			if(step == 0)
			fprintf(f,"%s  %s  %s  %s  %s  %s  %s  %s  %s  %s  %s  %s  %s  %s  %s  %s\n",
					"step_id","time","dt","circX","circY","circZ","linImpX","linImpY","linImpZ","angImpX","angImpY","angImpZ","Ens","Hel","Maxvor");

			fprintf(f, "%d  %9.9e  %9.9e  %9.9e  %9.9e  %9.9e  %9.9e  %9.9e  %9.9e  %9.9e  %9.9e  %9.9e  %9.9e  %9.9e  %9.9e\n",
					step, t, dt,globalSum[0],globalSum[1],globalSum[2],globalSum[3],globalSum[4],globalSum[5],
					globalSum[6],globalSum[7],globalSum[8],globalSum[10],globalSum[9],std::sqrt(maxVort));

			fclose(f);
		}
	}

	string getName() { return "Diagnostics"; }
};

#endif /* CoordinatorVorticity_h */
