//
//  CoordinatorPenalization.h
//  CubismUP_3D
//
//  Created by Christian Conti on 3/27/15.
//  Copyright (c) 2015 ETHZ. All rights reserved.
//

#ifndef CubismUP_3D_CoordinatorPenalization_h
#define CubismUP_3D_CoordinatorPenalization_h

#include "GenericOperator.h"
#include "GenericCoordinator.h"
//#include "OperatorPenalization.h"
//#include "Shape.h"

struct PenalizationObstacleVisitor : ObstacleVisitor
{
    FluidGrid * grid;
    const double dt, lambda, *uInf;
    vector<BlockInfo> vInfo;

    PenalizationObstacleVisitor(FluidGrid * grid, const double dt, const double lambda,const double* uInf)
    : grid(grid), dt(dt), lambda(lambda), uInf(uInf)
    {
        vInfo = grid->getBlocksInfo();
    }

     void visit(IF2D_ObstacleOperator* const obstacle)
     {
        obstacle->computeVelocities(uInf);

#pragma omp parallel
         {
            const std::map<int, ObstacleBlock*> obstblocks = obstacle->getObstacleBlocks();
            const Real factor = 0.5/(vInfo[0].h_gridpoint * dt);
            Real uBody[3], omegaBody[3], centerOfMass[3];
            obstacle->getCenterOfMass(centerOfMass);
            obstacle->getTranslationVelocity(uBody);
            obstacle->getAngularVelocity(omegaBody);

#pragma omp for schedule(static)
             for(int i=0; i<vInfo.size(); i++) {
                 const auto pos = obstblocks.find(i);
                 if(pos == obstblocks.end()) continue;

                 BlockInfo info = vInfo[i];
                 FluidBlock& b = *(FluidBlock*)info.ptrBlock;

				 for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
				 for(int iy=0; iy<FluidBlock::sizeY; ++iy)
				 for(int ix=0; ix<FluidBlock::sizeX; ++ix)
				 if (pos->second->chi[iz][iy][ix] > 0) {
					 Real p[3];
					 info.pos(p, ix, iy, iz);
					 p[0]-=centerOfMass[0];
					 p[1]-=centerOfMass[1];
					 p[2]-=centerOfMass[2];

					 const Real lamdtX  = dt * lambda * pos->second->chi[iz][iy][ix];
					 const Real object_UR[3] = {
							omegaBody[1]*p[2]-omegaBody[2]*p[1],
							omegaBody[2]*p[0]-omegaBody[0]*p[2],
							omegaBody[0]*p[1]-omegaBody[1]*p[0]
					 };
					 const Real object_UDEF[3] = {
							 pos->second->udef[iz][iy][ix][0],
							 pos->second->udef[iz][iy][ix][1],
							 pos->second->udef[iz][iy][ix][2]
					 };
					 b(ix,iy,iz).u = (b(ix,iy,iz).u + lamdtX * (uBody[0]+object_UR[0]+object_UDEF[0]-uInf[0])) / (1. + lamdtX);
					 b(ix,iy,iz).v = (b(ix,iy,iz).v + lamdtX * (uBody[1]+object_UR[1]+object_UDEF[1]-uInf[1])) / (1. + lamdtX);
					 b(ix,iy,iz).w = (b(ix,iy,iz).w + lamdtX * (uBody[2]+object_UR[2]+object_UDEF[2]-uInf[2])) / (1. + lamdtX);
             }
         }
     }
};

class CoordinatorPenalization : public GenericCoordinator
{
protected:
    IF2D_ObstacleVector** const obstacleVector;
    const double* const lambda;
    const double* const uInf;
public:
	CoordinatorPenalization(FluidGridMPI * grid, IF2D_ObstacleVector** const myobstacles, const double* const lambda, const double* const Uinf)
	: GenericCoordinator(grid), obstacleVector(myobstacles), lambda(lambda), uInf(Uinf)
	{ }
	
	void operator()(const double dt)
	{
		check((string)"penalization - start");

        PenalizationObstacleVisitor* penalizationVisitor = new PenalizationObstacleVisitor(grid, dt, *lambda, uInf);
        //std::vector<IF2D_ObstacleOperator*>* my_obstacles = (*obstacleVector)->getObstacleVector();
        //for(auto& obstacle : * my_obstacles) obstacle->Accept(penalizationVisitor);
        (*obstacleVector)->Accept(penalizationVisitor); // accept you son of a french cow
        delete penalizationVisitor;

		check((string)"penalization - end");
	}
	
	string getName()
	{
		return "Penalization";
	}
};

#endif
