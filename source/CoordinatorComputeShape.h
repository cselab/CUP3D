//
//  CubismUP_3D
//
//  Written by Guido Novati ( novatig@ethz.ch )
//  Copyright (c) 2017 ETHZ. All rights reserved.
//

#ifndef CubismUP_3D_CoordinatorComputeShape_h
#define CubismUP_3D_CoordinatorComputeShape_h

#include "GenericCoordinator.h"
#include "IF3D_ObstacleVector.h"

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
      for(int i=0; i<vInfo.size(); i++) {
        BlockInfo info = vInfo[i];
        FluidBlock& b = *(FluidBlock*)info.ptrBlock;

        for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
        for(int iy=0; iy<FluidBlock::sizeY; ++iy)
        for(int ix=0; ix<FluidBlock::sizeX; ++ix)
        b(ix,iy,iz).chi = 0;
      }
    }

    (*obstacleVector)->create(*stepID,*time, dt, Uinf);
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
  IF3D_ObstacleVector** obstacleVector;
  const double* const time;
  const Real* const Uinf;
  const double* const NU;
  const bool * const bDump;
  const int* const stepID;
public:
  CoordinatorComputeForces(FluidGridMPI*g, IF3D_ObstacleVector** myobst,
   const int* const s, const double* const _t, const double*const nu,
   const bool*const bdump, const Real*const uInf) : GenericCoordinator(g),
   obstacleVector(myobst), stepID(s), time(_t), NU(nu), bDump(bdump), Uinf(uInf)
    {    }

  void operator()(const double dt)
  {
    check((string)"obst. forces - start");
    (*obstacleVector)->computeForces(*stepID, *time, dt, Uinf, *NU, *bDump);
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
  IF3D_ObstacleVector** const obstacleVector;
  const double* const time;
  const double* const lambda;
  const Real* const Uinf;
  const int* const stepID;
public:
  CoordinatorComputeDiagnostics(FluidGridMPI*g, IF3D_ObstacleVector** myobst,
   const int* const s, const double* const _t, const double*const l,
   const Real*const u) : GenericCoordinator(g), obstacleVector(myobst), stepID(s), time(_t), lambda(l), Uinf(u) { }

  void operator()(const double dt)
  {
    check((string)"obst. diagnostics - start");
    (*obstacleVector)->computeDiagnostics(*stepID, *time, Uinf, *lambda);
    check((string)"obst. diagnostics - end");
  }

  string getName()
  {
    return "ComputeObstDiag";
  }
};

#endif
