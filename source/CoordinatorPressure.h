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
#include "GenericOperator.h"
#include "OperatorMovingPressure.h"
#include "PoissonSolverScalarFFTW_MPI.h"
#ifdef _MULTIGRID_
#include "MultigridHypre.h"
#endif // _MULTIGRID_

struct PressureObstacleVisitor : ObstacleVisitor
{
	FluidGridMPI * grid;
    vector<BlockInfo> vInfo;

    PressureObstacleVisitor(FluidGridMPI * grid):grid(grid)
    {
        vInfo = grid->getBlocksInfo();
    }

     void visit(const IF2D_ObstacleOperator* const obstacle)
     {
#pragma omp parallel
         {
             const std::map<int,ObstacleBlock*> obstblocks = obstacle->getObstacleBlocks();
             BlockInfo* const ary = &vInfo.front();
#pragma omp for schedule(static)
             for(int i=0; i<vInfo.size(); i++) {
                 const auto pos = obstblocks.find(i);
                 if(pos == obstblocks.end()) continue;

            	 FluidBlock& b = *(FluidBlock*)ary[i].ptrBlock;
                 for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
                 for(int iy=0; iy<FluidBlock::sizeY; ++iy)
				 for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
					 b(ix,iy,iz).tmpU += pos->second->udef[iz][iy][ix][0];
					 b(ix,iy,iz).tmpV += pos->second->udef[iz][iy][ix][1];
					 b(ix,iy,iz).tmpW += pos->second->udef[iz][iy][ix][2];
				 } //what if multiple obstacles share a block?? let's plus equal and worry about it
             }
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
		stencil = StencilInfo(-1,-1,-1, 2,2,2, false, 3, 1,2,3);
		
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
		for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
			const Real uW = lab(ix-1,iy  ,iz  ).u;
			const Real uE = lab(ix+1,iy  ,iz  ).u;
			const Real vS = lab(ix  ,iy-1,iz  ).v;
			const Real vN = lab(ix  ,iy+1,iz  ).v;
			const Real wF = lab(ix  ,iy  ,iz-1).w;
			const Real wB = lab(ix  ,iy  ,iz+1).w;
			o(ix, iy, iz).tmp  = factor * (uE-uW + vN-vS + wB-wF);
		}
	}
};

class OperatorDivergenceMinusDivTmpU : public GenericLabOperator
{
private:
	double dt;
	
public:
	OperatorDivergenceMinusDivTmpU(double dt) : dt(dt)
	{
		stencil = StencilInfo(-1,-1,-1, 2,2,2, false, 3, 1,2,3);
		
		stencil_start[0] = -1;
		stencil_start[1] = -1;
		stencil_start[2] = -1;
		stencil_end[0] = 2;
		stencil_end[1] = 2;
		stencil_end[2] = 2;
	}
	~OperatorDivergenceMinusDivTmpU() {}
	
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

			const Real uWdef = lab(ix-1,iy  ,iz  ).tmpU;
			const Real uEdef = lab(ix+1,iy  ,iz  ).tmpU;
			const Real vSdef = lab(ix  ,iy-1,iz  ).tmpV;
			const Real vNdef = lab(ix  ,iy+1,iz  ).tmpV;
			const Real wFdef = lab(ix  ,iy  ,iz-1).tmpW;
			const Real wBdef = lab(ix  ,iy  ,iz+1).tmpW;

			const Real divU = - factor * (uE - uW + vN - vS + wB - wF
					- lab(ix,iy).chi*(uEdef-uWdef+vNdef-vSdef+wBdef-wFdef));

			o(ix, iy, iz).tmp = divU;
			//o(ix, iy, iz).tmp  = divU;
		}
	}
};

class OperatorGradP : public GenericLabOperator
{
private:
	double dt;
	
public:
	OperatorGradP(double dt) : dt(dt)
	{
		stencil = StencilInfo(-1,-1,-1, 2,2,2, false, 2, 0,11);
		
		stencil_start[0] = -1;
		stencil_start[1] = -1;
		stencil_start[2] = -1;
		stencil_end[0] = 2;
		stencil_end[1] = 2;
		stencil_end[2] = 2;
	}
	
	~OperatorGradP() {}
	
	template <typename Lab, typename BlockType>
	void operator()(Lab & lab, const BlockInfo& info, BlockType& o) const
	{
		const double prefactor = -.5 * dt / (info.h_gridpoint);
		
		for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
		for(int iy=0; iy<FluidBlock::sizeY; ++iy)
		for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
			const Real divUW = lab(ix-1,iy  ,iz  ).tmp;
			const Real divUE = lab(ix+1,iy  ,iz  ).tmp;
			const Real divUS = lab(ix  ,iy-1,iz  ).tmp;
			const Real divUN = lab(ix  ,iy+1,iz  ).tmp;
			const Real divUF = lab(ix  ,iy  ,iz-1).tmp;
			const Real divUB = lab(ix  ,iy  ,iz+1).tmp;

			// divU contains the pressure correction after the Poisson solver
			o(ix,iy,iz).u += prefactor * (divUE - divUW);//
			o(ix,iy,iz).v += prefactor * (divUN - divUS);//
			o(ix,iy,iz).w += prefactor * (divUB - divUF);//

			assert(!std::isnan(o(ix,iy,iz).u));
			assert(!std::isnan(o(ix,iy,iz).v));
			assert(!std::isnan(o(ix,iy,iz).w));
		}
	}
};

template <typename Lab>
class CoordinatorPressure : public GenericCoordinator
{
protected:
    IF2D_ObstacleVector** const obstacleVector;
//#ifndef _MIXED_
	PoissonSolverScalarFFTW_MPI<FluidGridMPI, StreamerDiv> pressureSolver;
//#else
//	PoissonSolverScalarFFTW_MPI_DCT<FluidGridMPI, StreamerDiv> pressureSolver;
//#endif // _MIXED_
	
	inline void updatePressure()
	{
		const int N = vInfo.size();
		
#pragma omp parallel for schedule(static)
		for(int i=0; i<N; i++) {
			BlockInfo info = vInfo[i];
			FluidBlock& b = *(FluidBlock*)info.ptrBlock;
			
			for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
			for(int iy=0; iy<FluidBlock::sizeY; ++iy)
			for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
				b(ix,iy,iz).pOld = b(ix,iy,iz).p;
				b(ix,iy,iz).p    = b(ix,iy,iz).tmp;
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
			for (int i=0; i<N; i++) {
				lab.load(ary[i], 0);
				kernel(lab, ary[i], *(FluidBlock*)ary[i].ptrBlock);
			}
		}
	}
	
public:
	CoordinatorPressure(FluidGridMPI * grid, IF2D_ObstacleVector** const myobstacles) :
		GenericCoordinator(grid), pressureSolver(NTHREADS,*grid), obstacleVector(myobstacles)
	{
	}
	
	void operator()(const double dt)
	{
#pragma omp parallel for schedule(static)
		for(int i=0; i<vInfo.size(); i++) {
			FluidBlock& b = *(FluidBlock*)vInfo[i].ptrBlock;

			for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
			for(int iy=0; iy<FluidBlock::sizeY; ++iy)
			for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
				b(ix,iy).tmpU = 0;
				b(ix,iy).tmpV = 0;
				b(ix,iy).tmpW = 0; //zero fields, going to contain Udef
			}
		}

		PressureObstacleVisitor * pressureVisitor = new PressureObstacleVisitor(grid);
		//const auto my_obstacles = (*obstacleVector)->getObstacleVector();
		//for(const auto& obstacle : *my_obstacles) obstacle->Accept(pressureVisitor);
		(*obstacleVector)->Accept(penalizationVisitor); // accept you son of a french cow
		delete pressureVisitor;

		computeUnsplit<OperatorDivergenceMinusDivTmpU>(dt);
		pressureSolver.solve(*grid);
		computeUnsplit<OperatorGradP>(dt);
		updatePressure();
	}
	
	string getName()
	{
		return "Pressure";
	}
};
#endif
