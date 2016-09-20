//
//  CoordinatorPressure.h
//  CubismUP_3D
//
//  Created by Christian Conti on 3/30/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef CubismUP_3D_CoordinatorPressure_h
#define CubismUP_3D_CoordinatorPressure_h

#include "GenericOperator.h"
#include "GenericCoordinator.h"
#include "IF3D_ObstacleVector.h"
#include "PoissonSolverScalarFFTW_MPI.h"

struct PressureObstacleVisitor : ObstacleVisitor
{
	FluidGridMPI * grid;
    vector<BlockInfo> vInfo;

    PressureObstacleVisitor(FluidGridMPI * grid):grid(grid)
    {
        vInfo = grid->getBlocksInfo();
    }

     void visit(const IF3D_ObstacleOperator* const obstacle)
     {
#pragma omp parallel
         {
             const std::map<int,ObstacleBlock*> obstblocks = obstacle->getObstacleBlocks();
#pragma omp for schedule(static)
             for(int i=0; i<vInfo.size(); i++) {
            	 BlockInfo info = vInfo[i];
                 const auto pos = obstblocks.find(info.blockID);
                 if(pos == obstblocks.end()) continue;
                 FluidBlock& b = *(FluidBlock*)vInfo[i].ptrBlock;

                 for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
                 for(int iy=0; iy<FluidBlock::sizeY; ++iy)
				 for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
					 // what if multiple obstacles share a block??
					 // let's plus equal and wake up during the night to stress about it
					 b(ix,iy,iz).tmpU += pos->second->udef[iz][iy][ix][0];
					 b(ix,iy,iz).tmpV += pos->second->udef[iz][iy][ix][1];
					 b(ix,iy,iz).tmpW += pos->second->udef[iz][iy][ix][2];
				 }
             }
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
		stencil = StencilInfo(-1,-1,-1, 2,2,2, false, 6, 0,1,2,5,6,7);
		
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
		for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
			// Poisson solver reads field p for the rhs
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

			const Real divU = factor * (uE - uW + vN - vS + wB - wF
					-lab(ix,iy,iz).chi*(uEdef-uWdef+vNdef-vSdef+wBdef-wFdef) );
			o(ix, iy, iz).p = divU;
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
		stencil = StencilInfo(-1,-1,-1, 2,2,2, false, 1, 4);
		
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
			// p contains the pressure correction after the Poisson solver
			const Real divUW = lab(ix-1,iy  ,iz  ).p;
			const Real divUE = lab(ix+1,iy  ,iz  ).p;
			const Real divUS = lab(ix  ,iy-1,iz  ).p;
			const Real divUN = lab(ix  ,iy+1,iz  ).p;
			const Real divUF = lab(ix  ,iy  ,iz-1).p;
			const Real divUB = lab(ix  ,iy  ,iz+1).p;
			/*
			if(ix==0)
				printf("LB %d %d %d: %f | %f \n",info.blockID,iy,iz,lab(ix-1,iy  ,iz  ).p,lab(ix,iy  ,iz  ).p);
			if(ix==FluidBlock::sizeX-1)
				printf("UB %d %d %d: %f | %f \n",info.blockID,iy,iz,lab(ix,iy  ,iz  ).p,lab(ix+1,iy  ,iz  ).p);
			 */
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
    IF3D_ObstacleVector** const obstacleVector;
//#ifndef _MIXED_
	PoissonSolverScalarFFTW_MPI<FluidGridMPI, StreamerDiv> pressureSolver;
//#else
//	PoissonSolverScalarFFTW_MPI_DCT<FluidGridMPI, StreamerDiv> pressureSolver;
//#endif // _MIXED_

public:
	CoordinatorPressure(FluidGridMPI * grid, IF3D_ObstacleVector** const myobstacles) :
		GenericCoordinator(grid), pressureSolver(NTHREADS,*grid), obstacleVector(myobstacles)
	{
	}
	
	void operator()(const double dt)
	{
#pragma omp parallel
		{
			const int N = vInfo.size();

#pragma omp for schedule(static)
			for(int i=0; i<vInfo.size(); i++) {
				BlockInfo info = vInfo[i];
				FluidBlock& b = *(FluidBlock*)info.ptrBlock;
				for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
				for(int iy=0; iy<FluidBlock::sizeY; ++iy)
				for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
					b(ix,iy,iz).tmpU = 0;
					b(ix,iy,iz).tmpV = 0;
					b(ix,iy,iz).tmpW = 0; //zero fields, going to contain Udef
				}
			}
		}

		{	//store deformation velocities onto tmp fields
			ObstacleVisitor * pressureVisitor = new PressureObstacleVisitor(grid);
			(*obstacleVector)->Accept(pressureVisitor); // accept you son of a french cow
			delete pressureVisitor;
		}
		{ 	//place onto p: ( div u^(t+1) - div u^* ) / dt
			//where i want div u^(t+1) to be equal to div udef
			OperatorDivergenceMinusDivTmpU kernelDiv(dt);
			compute(kernelDiv);
		}
		{
			pressureSolver.solve(*grid);
		}
		{ //pressure correction dudt* = - grad P / rho
			OperatorGradP kernelGradP(dt);
			compute(kernelGradP);
		}
	}
	
	string getName()
	{
		return "Pressure";
	}
};
#endif
