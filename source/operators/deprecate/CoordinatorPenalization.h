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
#include "../obstacles/IF3D_ObstacleVector.h"

// for schooling: split first compute uinf (comp velocity with uinf = 0) and then penalize






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
