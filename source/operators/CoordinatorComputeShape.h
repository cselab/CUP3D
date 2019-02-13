//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#ifndef CubismUP_3D_CoordinatorComputeShape_h
#define CubismUP_3D_CoordinatorComputeShape_h

#include "GenericCoordinator.h"
#include "../obstacles/IF3D_ObstacleVector.h"

struct VelocityObstacleVisitor : public ObstacleVisitor
{
  FluidGridMPI * grid;
  const Real * const uInf;
  int * const nSum;
  double * const uSum;
  std::vector<BlockInfo> vInfo;
  bool finalize = false;

  VelocityObstacleVisitor(FluidGridMPI* _grid, const Real*const _uInf,
    int*const _nSum, double*const _uSum) : grid(_grid), uInf(_uInf), nSum(_nSum), uSum(_uSum)
  {
    vInfo = grid->getBlocksInfo();
  }

  void visit(IF3D_ObstacleOperator* const obstacle)
  {
    if (not finalize)
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
    else
    {
      double obstU[3];
      obstacle->getTranslationVelocity(obstU);
        obstU[0] += uInf[0]; //set obstacle speed to zero
        obstU[1] += uInf[1];
        obstU[2] += uInf[2];
      obstacle->setTranslationVelocity(obstU);
    }
  }
};

class CoordinatorComputeShape : public GenericCoordinator
{
 protected:
  IF3D_ObstacleVector** const obstacleVector;
  const double* const time;
  const Real* const Uinf;
  const int* const stepID;

 public:
  CoordinatorComputeShape(FluidGridMPI*g, IF3D_ObstacleVector**const myobst,
    const int*const s, const double*const _t, const Real*const uInf)
  : GenericCoordinator(g),obstacleVector(myobst),time(_t),Uinf(uInf),stepID(s)
  {
    //      (*obstacleVector)->create(*stepID,*time, 0, Uinf);
  }

  void operator()(const double dt)
  {
    check("shape - start");
    (*obstacleVector)->update(*stepID,*time, dt, Uinf);

    //each obstacle does chi = max(chi,obstacleblock->chi)
    #pragma omp parallel
    {
      #pragma omp for schedule(static)
      for (int i = 0; i < (int)vInfo.size(); ++i) {
        BlockInfo info = vInfo[i];
        FluidBlock& b = *(FluidBlock*)info.ptrBlock;

        for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
        for(int iy=0; iy<FluidBlock::sizeY; ++iy)
        for(int ix=0; ix<FluidBlock::sizeX; ++ix) b(ix,iy,iz).chi = 0;
      }
    }

    (*obstacleVector)->create(*stepID,*time, dt, Uinf);

    {
     int nSum[3] = {0,0,0};
     double uSum[3] = {0,0,0};
     ObstacleVisitor* velocityVisitor =
                     new VelocityObstacleVisitor(grid, uInf, nSum, uSum);
     (*obstacleVector)->Accept(velocityVisitor);//accept you son of a french cow
     if(nSum[0]) uInf[0] = uSum[0]/nSum[0];
     if(nSum[1]) uInf[1] = uSum[1]/nSum[1];
     if(nSum[2]) uInf[2] = uSum[2]/nSum[2];
     //printf("Old Uinf %g %g %g\n",uInf[0],uInf[1],uInf[2]);
     velocityVisitor->finalize = true;
     (*obstacleVector)->Accept(velocityVisitor);//accept you son of a french cow
     //if(rank == 0) if(nSum[0] || nSum[1] || nSum[2])
     //  printf("New Uinf %g %g %g (from %d %d %d)\n",
     //  uInf[0],uInf[1],uInf[2],nSum[0],nSum[1],nSum[2]);

     delete velocityVisitor;
    }
    check("shape - end");
  }

  std::string getName()
  {
    return "ComputeShape";
  }
};


class CoordinatorComputeForces : public GenericCoordinator
{
 protected:
  IF3D_ObstacleVector** obstacleVector;
  const int* const stepID;
  const double* const time;
  const double* const NU;
  const bool * const bDump;
  const Real* const Uinf;
public:
  CoordinatorComputeForces(FluidGridMPI*g, IF3D_ObstacleVector** myobst,
   const int* const s, const double* const _t, const double*const nu,
   const bool*const bdump, const Real*const uInf) : GenericCoordinator(g),
   obstacleVector(myobst), stepID(s), time(_t), NU(nu), bDump(bdump), Uinf(uInf)
    {    }

  void operator()(const double dt)
  {
    check((std::string)"obst. forces - start");
    (*obstacleVector)->computeForces(*stepID, *time, dt, Uinf, *NU, *bDump);
    check((std::string)"obst. forces - end");
  }

  std::string getName()
  {
      return "ComputeForces";
  }
};


class CoordinatorComputeDiagnostics : public GenericCoordinator
{
 protected:
  IF3D_ObstacleVector** const obstacleVector;
  const int* const stepID;
  const double* const time;
  const double* const lambda;
  const Real* const Uinf;
public:
  CoordinatorComputeDiagnostics(FluidGridMPI*g, IF3D_ObstacleVector** myobst,
   const int* const s, const double* const _t, const double*const l,
   const Real*const u) : GenericCoordinator(g), obstacleVector(myobst), stepID(s), time(_t), lambda(l), Uinf(u) { }

  void operator()(const double dt)
  {
    check((std::string)"obst. diagnostics - start");
    (*obstacleVector)->computeDiagnostics(*stepID, *time, Uinf, *lambda);
    check((std::string)"obst. diagnostics - end");
  }

  std::string getName()
  {
    return "ComputeObstDiag";
  }
};

#endif
