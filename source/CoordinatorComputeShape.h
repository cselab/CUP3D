//
//  CoordinatorComputeShape.h
//  CubismUP_3D
//
//  Created by Christian Conti on 3/30/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef CubismUP_3D_CoordinatorComputeShape_h
#define CubismUP_3D_CoordinatorComputeShape_h

#include "GenericCoordinator.h"
#include "OperatorComputeShape.h"
#include "Shape.h"

class CoordinatorComputeShape : public GenericCoordinator
{
protected:
	IF2D_ObstacleVector** const obstacleVector;
    const double* const time;
    const double* const Uinf;
    const int* const stepID;
    
public:
    CoordinatorComputeShape(FluidGrid * grid, IF2D_ObstacleVector** const myobstacles, const int* const stepID, const double * const time, const double * const Uinf)
    : GenericCoordinator(grid), obstacleVector(myobstacles), time(time), Uinf(Uinf), stepID(stepID)
	{
    	//std::vector<IF2D_ObstacleOperator*>* my_obstacles = (*obstacleVector)->getObstacleVector();
    	//for(const auto & obstacle_ptr : * my_obstacles) obstacle_ptr->create(*stepID, *time, 0);
    	(*obstacleVector)->create(*stepID,*time, 0);
	}
	
	void operator()(const double dt)
	{
		check("shape - start");
		
        //std::vector<IF2D_ObstacleOperator*>* my_obstacles = (*obstacleVector)->getObstacleVector();
        //for(const auto & obstacle_ptr : * my_obstacles) obstacle_ptr->update(*stepID,*time, dt);
		(*obstacleVector)->update(*stepID,*time, dt);

#pragma omp parallel
        {
#pragma omp for schedule(static)
        	for(int i=0; i<vInfo.size(); i++) {
        		BlockInfo info = vInfo[i];
        		FluidBlock& b = *(FluidBlock*)info.ptrBlock;

                for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
        		for(int iy=0; iy<FluidBlock::sizeY; ++iy)
				for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
					b(ix,iy,iz).chi = 0;
					b(ix,iy,iz).tmp = 0;
				}
        	}
        }

        //for(const auto & obstacle_ptr : * my_obstacles) obstacle_ptr->create(*stepID, *time, dt);
        (*obstacleVector)->create(*stepID,*time, dt);
		check("shape - end");
	}
	
	string getName()
	{
		return "ComputeShape";
	}
};


class CoordinatorComputeForces : public GenericCoordinator
{
protected:
    IF2D_ObstacleVector** obstacleVector;
    const double* const time;
    const double* const Uinf;
    const double* const NU;
    const bool * const bDump;
    const int* const stepID;
public:
    CoordinatorComputeForces(FluidGrid * grid, IF2D_ObstacleVector** myobstacles, const int* const stepID, const double* const time, const double* const NU, const bool* const bDump, const double* const Uinf)
: GenericCoordinator(grid), obstacleVector(myobstacles), stepID(stepID), time(time), NU(NU), bDump(bDump), Uinf(Uinf)
    {    }

    void operator()(const double dt)
    {
        check((string)"obst. forces - start");

        //std::vector<IF2D_ObstacleOperator*>* my_obstacles = (*obstacleVector)->getObstacleVector();
        //for(const auto & obstacle_ptr : * my_obstacles) obstacle_ptr->computeForces(*stepID, *time, Uinf, *NU, *bDump);
        (*obstacleVector)->computeForces(*stepID, *time, Uinf, *NU, *bDump);
        check((string)"obst. forces - end");
    }

    string getName()
    {
        return "ComputeForces";
    }
};


class CoordinatorComputeDiagnostics : public GenericCoordinator
{
protected:
    IF2D_ObstacleVector** const obstacleVector;
    const double* const time;
    const double* const lambda;
    const double* const Uinf;
    const int* const stepID;
public:
    CoordinatorComputeDiagnostics(FluidGrid * grid, IF2D_ObstacleVector** const myobstacles, const int* const stepID, const double* const time, const double* const lambda, const double* const Uinf)
: GenericCoordinator(grid), obstacleVector(myobstacles), stepID(stepID), time(time), lambda(lambda), Uinf(Uinf)
    {    }

    void operator()(const double dt)
    {
        check((string)"obst. diagnostics - start");

        //std::vector<IF2D_ObstacleOperator*>* my_obstacles = (*obstacleVector)->getObstacleVector();
        //for(const auto & obstacle_ptr : * my_obstacles) obstacle_ptr->computeDiagnostics(*stepID, *time, Uinf, *lambda);
        (*obstacleVector)->computeDiagnostics(*stepID, *time, Uinf, *bDump);
        check((string)"obst. diagnostics - end");
    }

    string getName()
    {
        return "ComputeObstDiag";
    }
};

#endif
