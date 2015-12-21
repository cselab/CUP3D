//
//  CoordinatorPressure.h
//  CubismUP_3D
//
//  Created by Christian Conti on 3/30/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef CubismUP_3D_CoordinatorPressure_h
#define CubismUP_3D_CoordinatorPressure_h

#include "GenericCoordinator.h"
#include "OperatorDivergence.h"
#include "OperatorGradP.h"
#include "PoissonSolverScalarFFTW.h"
#ifdef _MULTIGRID_
#include "MultigridHypre.h"
#endif // _MULTIGRID_

#define _HYDROSTATIC_

template <typename Lab>
class CoordinatorPressure : public GenericCoordinator
{
protected:
	const int rank, nprocs;
	const double minRho;
	Real gravity[3];
	int * step;
	const bool bSplit;
	Real *uBody, *vBody, *wBody;
	
#ifdef _SPLIT_
#ifndef _MIXED_
	PoissonSolverScalarFFTW<FluidGrid, StreamerDiv> pressureSolver;
#else // _MIXED_
	PoissonSolverScalarFFTW_DCT<FluidGrid, StreamerDiv> pressureSolver;
#endif // _MIXED_
#endif // _SPLIT_
	
#ifdef _MULTIGRID_
	MultigridHypre mg;
#endif // _MULTIGRID_
	
	inline void addHydrostaticPressure(const double dt)
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
						b(ix,iy,iz).u -= dt*gravity[0]/b(ix,iy,iz).rho;
						b(ix,iy,iz).v -= dt*gravity[1]/b(ix,iy,iz).rho;
						b(ix,iy,iz).w -= dt*gravity[2]/b(ix,iy,iz).rho;
					}
		}
	}
	
	inline void updatePressure()
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
						b(ix,iy,iz).pOld = b(ix,iy,iz).p;
						b(ix,iy,iz).p    = b(ix,iy,iz).divU;
					}
		}
	}
	
	template <typename Operator>
	void computeSplit(const double dt)
	{
		Operator kernel(dt, minRho, *step);
		compute(kernel);
		/*
		BlockInfo * ary = &vInfo.front();
		const int N = vInfo.size();
		
#pragma omp parallel
		{
			Operator kernel(dt, minRho, *step);
			
			Lab lab;
			lab.prepare(*grid, kernel.stencil_start, kernel.stencil_end, false);
			
#pragma omp for schedule(static)
			for (int i=0; i<N; i++)
			{
				lab.load(ary[i], 0);
				kernel(lab, ary[i], *(FluidBlock*)ary[i].ptrBlock);
			}
		}
		*/
	}
	
#ifdef _SPLIT_
	template <typename Operator>
	void computeSplitFFTW(const double dt)
	{
		BlockInfo * ary = &vInfo.front();
		const int N = vInfo.size();
		
#pragma omp parallel
		{
			int dim[3] = {grid->getBlocksPerDimension(0)*FluidBlock::sizeX, grid->getBlocksPerDimension(1)*FluidBlock::sizeY, grid->getBlocksPerDimension(2)*FluidBlock::sizeZ};
			Operator kernel(dt, minRho, *step, pressureSolver.data, dim);
			
			Lab lab;
			lab.prepare(*grid, kernel.stencil_start, kernel.stencil_end, false);
			
#pragma omp for schedule(static)
			for (int i=0; i<N; i++)
			{
				lab.load(ary[i], 0);
				kernel(lab, ary[i], *(FluidBlock*)ary[i].ptrBlock);
			}
		}
	}
#endif
	
	template <typename Operator>
	void computeUnsplit(const double dt)
	{
		Operator kernel(dt);
		compute(kernel);
		/*
		BlockInfo * ary = &vInfo.front();
		const int N = vInfo.size();
		
#pragma omp parallel
		{
			Operator kernel(dt);
			
			Lab lab;
			lab.prepare(*grid, kernel.stencil_start, kernel.stencil_end, true);
			
#pragma omp for schedule(static)
			for (int i=0; i<N; i++)
			{
				lab.load(ary[i], 0);
				
				kernel(lab, ary[i], *(FluidBlock*)ary[i].ptrBlock);
			}
		}
		 */
	}
	
public:
	CoordinatorPressure(const double minRho, const Real gravity[3], Real * uBody, Real * vBody, Real * wBody, int * step, const bool bSplit, FluidGrid * grid, const int rank, const int nprocs) : GenericCoordinator(grid), rank(rank), nprocs(nprocs), minRho(minRho), step(step), bSplit(bSplit), uBody(uBody), vBody(vBody), wBody(wBody), gravity{gravity[0],gravity[1],gravity[2]}
#ifdef _SPLIT_
	, pressureSolver(NTHREADS,*grid)
#endif // _SPLIT_
	{
	}
	
	CoordinatorPressure(const double minRho, const Real gravity[3], int * step, const bool bSplit, FluidGrid * grid, const int rank, const int nprocs) : GenericCoordinator(grid), rank(rank), nprocs(nprocs), minRho(minRho), step(step), bSplit(bSplit), uBody(NULL), vBody(NULL), gravity{gravity[0],gravity[1],gravity[2]}
#ifdef _SPLIT_
	, pressureSolver(NTHREADS,*grid)
#endif // _SPLIT_
	{
	}
	
	void operator()(const double dt)
	{
		// need an interface that is the same for all solvers - this way the defines can be removed more cleanly
		MPI_Barrier(MPI_COMM_WORLD);
		
		check("pressure - start");
		
		// pressure
#ifdef _SPLIT_
#ifdef _HYDROSTATIC_
		addHydrostaticPressure(dt);
#endif // _HYDROSTATIC_
		computeSplit<OperatorDivergenceSplit>(dt); // this part could be done directly in the correct data structure
		//computeSplitFFTW<OperatorDivergenceSplitFFTW>(dt); // this part could be done directly in the correct data structure
		pressureSolver.solve(*grid,false);
		computeSplit<OperatorGradPSplit>(dt); // this part could be done directly in the correct data structure
#endif // _SPLIT_
#ifdef _MULTIGRID_
		if (rank==0)
		{
#ifdef _HYDROSTATIC_
			addHydrostaticPressure(dt);
#endif
			if (bSplit)
				computeSplit<OperatorDivergenceSplit>(dt);
			else
				compute<OperatorDivergence>(dt);
		}
		check("pressure - preMG");
		mg.setup(grid, bSplit, rank, nprocs);
		mg();
		
		MPI_Barrier(MPI_COMM_WORLD);
		
		check("pressure - postMG");
		if (rank==0)
		{
			if (bSplit)
				computeSplit<OperatorGradPSplit>(dt);
			else
				computeUnsplit<OperatorGradP>(dt);
		}
#endif // _MULTIGRID_
		
		updatePressure();
		
		check("pressure - end");
	}
	
	string getName()
	{
		return "Pressure";
	}
};

template <typename Lab>
class CoordinatorPressureSimple : public GenericCoordinator
{
protected:
#ifndef _MIXED_
	PoissonSolverScalarFFTW<FluidGrid, StreamerDiv> pressureSolver;
#else
	PoissonSolverScalarFFTW_DCT<FluidGrid, StreamerDiv> pressureSolver;
#endif // _MIXED_
	
	inline void updatePressure()
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
						b(ix,iy,iz).pOld = b(ix,iy,iz).p;
						b(ix,iy,iz).p    = b(ix,iy,iz).divU;
					}
		}
	}
	
	template <typename Operator>
	void computeUnsplit(const double dt)
	{
		BlockInfo * ary = &vInfo.front();
		const int N = vInfo.size();
		
#pragma omp parallel
		{
			Operator kernel(dt);
			
			Lab lab;
			lab.prepare(*grid, kernel.stencil_start, kernel.stencil_end, true);
			
#pragma omp for schedule(static)
			for (int i=0; i<N; i++)
			{
				lab.load(ary[i], 0);
				kernel(lab, ary[i], *(FluidBlock*)ary[i].ptrBlock);
			}
		}
	}
	
public:
	CoordinatorPressureSimple(FluidGrid * grid) : GenericCoordinator(grid), pressureSolver(NTHREADS,*grid)
	{
	}
	
	void operator()(const double dt)
	{
		computeUnsplit<OperatorDivergence>(dt);
		pressureSolver.solve(*grid,true);
		computeUnsplit<OperatorGradP>(dt);
		
		updatePressure();
	}
	
	string getName()
	{
		return "Pressure";
	}
};
#endif
