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
	Real gravity[2];
	int * step;
    const bool bSplit;
    Real *uBody, *vBody;
	
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
			
			for(int iy=0; iy<FluidBlock::sizeY; ++iy)
				for(int ix=0; ix<FluidBlock::sizeX; ++ix)
				{
					b(ix,iy).u -= dt*gravity[0]/b(ix,iy).rho;
					b(ix,iy).v -= dt*gravity[1]/b(ix,iy).rho;
					
					// doesn't help much
					//b(ix,iy).u -= dt*minRho*gravity[0]/b(ix,iy).rho + (minRho<1 ? dt*(1-minRho)*gravity[0] : 0);
					//b(ix,iy).v -= dt*minRho*gravity[1]/b(ix,iy).rho + (minRho<1 ? dt*(1-minRho)*gravity[1] : 0);
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
			
			for(int iy=0; iy<FluidBlock::sizeY; ++iy)
				for(int ix=0; ix<FluidBlock::sizeX; ++ix)
				{
					b(ix,iy).pOld = b(ix,iy).p;
					b(ix,iy).p    = b(ix,iy).divU;
				}
		}
	}
	
	template <typename Operator>
	void computeSplit(const double dt)
	{
		BlockInfo * ary = &vInfo.front();
		const int N = vInfo.size();
		
#pragma omp parallel
		{
			Operator kernel(dt, minRho, *step);
			
			Lab mylab;
			mylab.prepare(*grid, kernel.stencil_start, kernel.stencil_end, false);
			
#pragma omp for schedule(static)
			for (int i=0; i<N; i++)
			{
				mylab.load(ary[i], 0);
				kernel(mylab, ary[i], *(FluidBlock*)ary[i].ptrBlock);
			}
		}
	}
	
	template <typename Operator>
	void compute(const double dt)
	{
		BlockInfo * ary = &vInfo.front();
		const int N = vInfo.size();
		
#pragma omp parallel
		{
			Operator kernel(dt);
			
            Lab mylab;
#ifdef _MOVING_FRAME_
            mylab.pDirichlet.u = 0;
            mylab.pDirichlet.v = *vBody;
#endif
			mylab.prepare(*grid, kernel.stencil_start, kernel.stencil_end, true);
			
#pragma omp for schedule(static)
			for (int i=0; i<N; i++)
			{
				mylab.load(ary[i], 0);
				
				kernel(mylab, ary[i], *(FluidBlock*)ary[i].ptrBlock);
			}
		}
	}
	
public:
	CoordinatorPressure(const double minRho, const Real gravity[2], Real * uBody, Real * vBody, int * step, const bool bSplit, FluidGrid * grid, const int rank, const int nprocs) : GenericCoordinator(grid), rank(rank), nprocs(nprocs), minRho(minRho), step(step), bSplit(bSplit), uBody(uBody), vBody(vBody), gravity{gravity[0],gravity[1]}
#ifdef _SPLIT_
	, pressureSolver(NTHREADS,*grid)
#endif // _SPLIT_
	{
	}
    
    CoordinatorPressure(const double minRho, const Real gravity[2], int * step, const bool bSplit, FluidGrid * grid, const int rank, const int nprocs) : GenericCoordinator(grid), rank(rank), nprocs(nprocs), minRho(minRho), step(step), bSplit(bSplit), uBody(NULL), vBody(NULL), gravity{gravity[0],gravity[1]}
#ifdef _SPLIT_
    , pressureSolver(NTHREADS,*grid)
#endif // _SPLIT_
    {
    }
	
	void operator()(const double dt)
	{
		// need an interface that is the same for all solvers - this way the defines can be removed more cleanly
#ifdef _MULTIGRID_
		MPI_Barrier(MPI_COMM_WORLD);
#endif // _MULTIGRID_
		
		check("pressure - start");
		
		// pressure
#ifdef _SPLIT_
#ifdef _HYDROSTATIC_
		addHydrostaticPressure(dt);
#endif // _HYDROSTATIC_
		computeSplit<OperatorDivergenceSplit>(dt);
		pressureSolver.solve(*grid,false);
		computeSplit<OperatorGradPSplit>(dt);
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
				compute<OperatorGradP>(dt);
		}
#endif // _MULTIGRID_
		
#ifdef _MULTIGRID_
		if (rank==0)
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
			
			for(int iy=0; iy<FluidBlock::sizeY; ++iy)
				for(int ix=0; ix<FluidBlock::sizeX; ++ix)
				{
					b(ix,iy).pOld = b(ix,iy).p;
					b(ix,iy).p    = b(ix,iy).divU;
				}
		}
	}
	
	template <typename Operator>
	void compute(const double dt)
	{
		BlockInfo * ary = &vInfo.front();
		const int N = vInfo.size();
		
#pragma omp parallel
		{
			Operator kernel(dt);
			
			Lab mylab;
			mylab.prepare(*grid, kernel.stencil_start, kernel.stencil_end, true);
			
#pragma omp for schedule(static)
			for (int i=0; i<N; i++)
			{
				mylab.load(ary[i], 0);
				
				kernel(mylab, ary[i], *(FluidBlock*)ary[i].ptrBlock);
			}
		}
	}
	
public:
	CoordinatorPressureSimple(FluidGrid * grid) : GenericCoordinator(grid), pressureSolver(NTHREADS,*grid)
	{
	}
	
	void operator()(const double dt)
	{
		compute<OperatorDivergence>(dt);
		pressureSolver.solve(*grid,true);
		compute<OperatorGradP>(dt);
		
		updatePressure();
	}
	
	string getName()
	{
		return "Pressure";
	}
};
#endif
