//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
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
  const double dt, lambda;
  const Real * const uInf;
  //Real ext_X, ext_Y, ext_Z;
  std::vector<BlockInfo> vInfo;

  PenalizationObstacleVisitor(FluidGridMPI* g, const double _dt,
      const double l, const Real*const u) : grid(g), dt(_dt),
    lambda(l), uInf(u)
  {
    vInfo = grid->getBlocksInfo();
  }

  void visit(IF3D_ObstacleOperator* const obstacle)
  {
    {
      double obstU[3];
      obstacle->getTranslationVelocity(obstU);
      //if (obstacle->obstacleID<2) {
      //  if(!obstacle->rank)
      //    printf("Discrepancy of obstacle %d = %g %g\n",
      //            obstacle->obstacleID, obstU[0]+uInf[0], obstU[1]+uInf[1]);
      //  obstU[0] = 0;
      //  obstU[1] = 0;
      //} else {
        obstU[0] += uInf[0]; //set obstacle speed to zero
        obstU[1] += uInf[1];
      //}
        obstU[2] += uInf[2];
      obstacle->setTranslationVelocity(obstU);
    }

    #pragma omp parallel
    {
      const std::map<int, ObstacleBlock*>& obstblocks = obstacle->getObstacleBlocks();
      double uBody[3], omegaBody[3], centerOfMass[3];
      obstacle->getCenterOfMass(centerOfMass);
      obstacle->getTranslationVelocity(uBody);
      obstacle->getAngularVelocity(omegaBody);
      const Real lamdt = double(dt) * double(lambda);
      #pragma omp for schedule(dynamic)
      for (int i = 0; i < (int)vInfo.size(); ++i) {
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
          const Real lamdtX = lamdt * pos->second->chi[iz][iy][ix];
          const Real object_UR[3] = {
              (Real)omegaBody[1]*p[2] -(Real)omegaBody[2]*p[1],
              (Real)omegaBody[2]*p[0] -(Real)omegaBody[0]*p[2],
              (Real)omegaBody[0]*p[1] -(Real)omegaBody[1]*p[0]
          };
          const Real object_UDEF[3] = {
              pos->second->udef[iz][iy][ix][0],
              pos->second->udef[iz][iy][ix][1],
              pos->second->udef[iz][iy][ix][2]
          };
          const Real U_TOT[3] = {
              (Real)uBody[0] +object_UR[0] +object_UDEF[0] -uInf[0],
              (Real)uBody[1] +object_UR[1] +object_UDEF[1] -uInf[1],
              (Real)uBody[2] +object_UR[2] +object_UDEF[2] -uInf[2]
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
  const Real * const uInf;
  int * const nSum;
  double * const uSum;
  std::vector<BlockInfo> vInfo;

  VelocityObstacleVisitor(FluidGridMPI* _grid, const Real*const _uInf,
    int*const _nSum, double*const _uSum) : grid(_grid), uInf(_uInf), nSum(_nSum), uSum(_uSum)
  {
    vInfo = grid->getBlocksInfo();
  }

  void visit(IF3D_ObstacleOperator* const obstacle)
  {
    const auto &bFixFrameOfRef = obstacle->bFixFrameOfRef;
    const Real dummy[3] = { 0.0, 0.0, 0.0 };
    obstacle->computeVelocities(dummy); // compute velocities with zero uinf
    double povU[3];
    obstacle->getTranslationVelocity(povU);

    if (bFixFrameOfRef[0]) { (nSum[0])++; uSum[0] -= povU[0]; }
    if (bFixFrameOfRef[1]) { (nSum[1])++; uSum[1] -= povU[1]; }
    if (bFixFrameOfRef[2]) { (nSum[2])++; uSum[2] -= povU[2]; }
  }
};

class CoordinatorPenalization : public GenericCoordinator
{
protected:
    IF3D_ObstacleVector** const obstacleVector;
    double* const lambda;
    Real* const uInf;
    int rank = 0;
public:
  CoordinatorPenalization(FluidGridMPI * g, IF3D_ObstacleVector** const myobst, double* const l, Real* const u)
  : GenericCoordinator(g), obstacleVector(myobst), lambda(l), uInf(u)
  {
    MPI_Comm comm = g->getCartComm();
    MPI_Comm_rank(comm, &rank);
  }

  void operator()(const double dt)
  {
    check((std::string)"penalization - start");

    int nSum[3] = {0,0,0};
    double uSum[3] = {0,0,0};
    ObstacleVisitor* velocityVisitor =
                    new VelocityObstacleVisitor(grid, uInf, nSum, uSum);
    (*obstacleVector)->Accept(velocityVisitor); // accept you son of a french cow
    if(nSum[0]) uInf[0] = uSum[0]/nSum[0];
    if(nSum[1]) uInf[1] = uSum[1]/nSum[1];
    if(nSum[2]) uInf[2] = uSum[2]/nSum[2];
      //printf("Old Uinf %g %g %g\n",uInf[0],uInf[1],uInf[2]);

    //if(rank == 0) if(nSum[0] || nSum[1] || nSum[2])
    //  printf("New Uinf %g %g %g (from %d %d %d)\n",
    //  uInf[0],uInf[1],uInf[2],nSum[0],nSum[1],nSum[2]);

    delete velocityVisitor;

    ObstacleVisitor* penalizationVisitor =
                    new PenalizationObstacleVisitor(grid, dt, *lambda, uInf);
    (*obstacleVector)->Accept(penalizationVisitor); // accept you son of a french cow
    delete penalizationVisitor;

    check((std::string)"penalization - end");
  }

  std::string getName()
  {
    return "Penalization";
  }
};

#endif
