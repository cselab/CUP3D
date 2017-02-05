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
#include "IF3D_ObstacleVector.h"

// for schooling: split first compute uinf (comp velocity with uinf = 0) and then penalize
struct PenalizationObstacleVisitor : public ObstacleVisitor
{
	FluidGridMPI * grid;
  const Real dt, lambda;
  Real * const uInf;
  //Real ext_X, ext_Y, ext_Z;
  vector<BlockInfo> vInfo;

  PenalizationObstacleVisitor(FluidGridMPI * grid, const Real dt,
															const Real lambda, Real* const uInf)
  : grid(grid), dt(dt), lambda(lambda), uInf(uInf)
  {
    vInfo = grid->getBlocksInfo();
  }

 void visit(IF3D_ObstacleOperator* const obstacle)
 {
	 const bool bFixFrameOfRef = obstacle->bFixFrameOfRef;
	 if (bFixFrameOfRef) {
		 if (obstacle->obstacleID!=0) {
			 printf("Can only fix first obstacle.\n"); abort();
		 }
		 // compute velocities with zero uinf
     Real leadU[3], dummy[3] = {0.0, 0.0, 0.0};
     obstacle->computeVelocities(dummy);
     obstacle->getTranslationVelocity(leadU);
     uInf[0] = -leadU[0]; //uInf now is speed of this obstacle
     uInf[1] = -leadU[1];
     uInf[2] = -leadU[2];
     leadU[0] += uInf[0]; //set obstacle speed to zero
     leadU[1] += uInf[1];
     leadU[2] += uInf[2];
     obstacle->setTranslationVelocity(leadU);
	 } else {
		 obstacle->computeVelocities(uInf);
	 }

	 #pragma omp parallel
   {
      const std::map<int, ObstacleBlock*> obstblocks = obstacle->getObstacleBlocks();
      Real uBody[3], omegaBody[3], centerOfMass[3];
      obstacle->getCenterOfMass(centerOfMass);
      obstacle->getTranslationVelocity(uBody);
      obstacle->getAngularVelocity(omegaBody);
      /*
			#pragma omp master
      printf("%f %f %f %f %f %f %f %f %f %f %f %f\n",
      		uBody[0],uBody[1],uBody[2],omegaBody[0],omegaBody[1],omegaBody[2],
      		centerOfMass[0],centerOfMass[1],centerOfMass[2],uInf[0],uInf[1],uInf[2]);
					*/
			#pragma omp for schedule(dynamic)
      for(int i=0; i<vInfo.size(); i++) {
      	BlockInfo info = vInfo[i];
      	const auto pos = obstblocks.find(info.blockID);
      	if(pos == obstblocks.end()) continue;
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
					const Real U_TOT[3] = {
							uBody[0]+object_UR[0]+object_UDEF[0]-uInf[0],
							uBody[1]+object_UR[1]+object_UDEF[1]-uInf[1],
							uBody[2]+object_UR[2]+object_UDEF[2]-uInf[2]
					};
					const Real alpha = 1./(1.+lamdtX);
					b(ix,iy,iz).u = alpha*b(ix,iy,iz).u + (1.-alpha)*(U_TOT[0]);
					b(ix,iy,iz).v = alpha*b(ix,iy,iz).v + (1.-alpha)*(U_TOT[1]);
					b(ix,iy,iz).w = alpha*b(ix,iy,iz).w + (1.-alpha)*(U_TOT[2]);
				}
     }
   }
 }
};

struct VelocityObstacleVisitor : public ObstacleVisitor
{
	FluidGridMPI * grid;
	Real*const uSum;
	int*const nSum;
  const Real * const uInf;
  vector<BlockInfo> vInfo;

    VelocityObstacleVisitor(FluidGridMPI * _grid, const Real* const _uInf,
														int*const _nSum, Real*const _uSum)
    : grid(_grid), uInf(_uInf), nSum(_nSum), uSum(_uSum)
    {
        vInfo = grid->getBlocksInfo();
    }

   void visit(IF3D_ObstacleOperator* const obstacle)
   {
		 const bool bFixFrameOfRef = obstacle->bFixFrameOfRef;
		 const Real dummy[3] = { 0.0, 0.0, 0.0 };
		 obstacle->computeVelocities(dummy);
  	 if (bFixFrameOfRef) {
			 (*nSum)++;
			 // compute velocities with zero uinf
       Real povU[3];
       obstacle->getTranslationVelocity(povU);
       uSum[0] -= povU[0]; //uInf now is speed of this obstacle
       uSum[1] -= povU[1];
       uSum[2] -= povU[2];
  	 }
	 }
};

struct POVPenalizationObstacleVisitor : public ObstacleVisitor
{
	FluidGridMPI * grid;
  const Real dt, lambda;
  const Real * const uInf;
  //Real ext_X, ext_Y, ext_Z;
  vector<BlockInfo> vInfo;

  POVPenalizationObstacleVisitor(FluidGridMPI * grid, const Real dt,
																const Real lambda, const Real* const uInf)
  : grid(grid), dt(dt), lambda(lambda), uInf(uInf)
  {
    vInfo = grid->getBlocksInfo();
  }

 void visit(IF3D_ObstacleOperator* const obstacle)
 {
	 // compute velocities with zero uinf
   Real obstU[3];
   obstacle->getTranslationVelocity(obstU);
	if (obstacle->obstacleID<2) {
		if(!obstacle->rank)
			printf("Discrepancy of obstacle %d = %g %g\n",
							obstacle->obstacleID, obstU[0]+uInf[0], obstU[1]+uInf[1]);
		obstU[0] = 0;
		obstU[1] = 0;
	} else {
	  obstU[0] += uInf[0]; //set obstacle speed to zero
	  obstU[1] += uInf[1];
	}
   obstU[2] += uInf[2];
   obstacle->setTranslationVelocity(obstU);

	 #pragma omp parallel
   {
      const std::map<int, ObstacleBlock*> obstblocks = obstacle->getObstacleBlocks();
      Real uBody[3], omegaBody[3], centerOfMass[3];
      obstacle->getCenterOfMass(centerOfMass);
      obstacle->getTranslationVelocity(uBody);
      obstacle->getAngularVelocity(omegaBody);

			#pragma omp for schedule(dynamic)
      for(int i=0; i<vInfo.size(); i++) {
      	BlockInfo info = vInfo[i];
      	const auto pos = obstblocks.find(info.blockID);
      	if(pos == obstblocks.end()) continue;
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
					const Real U_TOT[3] = {
							uBody[0]+object_UR[0]+object_UDEF[0]-uInf[0],
							uBody[1]+object_UR[1]+object_UDEF[1]-uInf[1],
							uBody[2]+object_UR[2]+object_UDEF[2]-uInf[2]
					};
					const Real alpha = 1./(1.+lamdtX);
					b(ix,iy,iz).u = alpha*b(ix,iy,iz).u + (1.-alpha)*(U_TOT[0]);
					b(ix,iy,iz).v = alpha*b(ix,iy,iz).v + (1.-alpha)*(U_TOT[1]);
					b(ix,iy,iz).w = alpha*b(ix,iy,iz).w + (1.-alpha)*(U_TOT[2]);
				}
     }
   }
 }
};

class CoordinatorPenalization : public GenericCoordinator
{
protected:
    IF3D_ObstacleVector** const obstacleVector;
    Real* const lambda;
    Real* const uInf;
public:
	CoordinatorPenalization(FluidGridMPI * grid, IF3D_ObstacleVector** const myobstacles, Real* const lambda, Real* const Uinf)
	: GenericCoordinator(grid), obstacleVector(myobstacles), lambda(lambda), uInf(Uinf)
	{ }

	void operator()(const Real dt)
	{
		check((string)"penalization - start");
#ifdef __2Leads_
		int nSum = 0;
		Real uSum[3] = {0,0,0};
    ObstacleVisitor* velocityVisitor =
										new VelocityObstacleVisitor(grid, uInf, &nSum, uSum);
    (*obstacleVector)->Accept(velocityVisitor); // accept you son of a french cow
		if(nSum) {
			printf("Old Uinf %g %g %g\n",uInf[0],uInf[1],uInf[2]);
			uInf[0] = uSum[0]/nSum;
			uInf[1] = uSum[1]/nSum;
			uInf[2] = uSum[2]/nSum;
			printf("New Uinf %g %g %g (from %d)\n",uInf[0],uInf[1],uInf[2],nSum);
		}
    delete velocityVisitor;

    ObstacleVisitor* penalizationVisitor =
										new POVPenalizationObstacleVisitor(grid, dt, *lambda, uInf);
    (*obstacleVector)->Accept(penalizationVisitor); // accept you son of a french cow
    delete penalizationVisitor;
#else
    ObstacleVisitor* penalizationVisitor =
											new PenalizationObstacleVisitor(grid, dt, *lambda, uInf);
    (*obstacleVector)->Accept(penalizationVisitor); // accept you son of a french cow
    delete penalizationVisitor;
#endif
		check((string)"penalization - end");
	}

	string getName()
	{
		return "Penalization";
	}
};

#endif
