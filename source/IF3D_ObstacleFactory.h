//
//  CubismUP_3D
//
//  Written by Guido Novati ( novatig@ethz.ch ).
//  This file started as an extension of code written by Wim van Rees
//  Copyright (c) 2017 ETHZ. All rights reserved.
//

#ifndef __IF3D_ROCKS__IF3D_ObstacleFactory__
#define __IF3D_ROCKS__IF3D_ObstacleFactory__

#include "IF3D_ObstacleOperator.h"

class IF3D_ObstacleFactory
{
  int rank;
  FluidGridMPI * grid;
  const Real * const Uinf;
  int _getlines(std::string filename);

public:
  IF3D_ObstacleFactory(FluidGridMPI*g, const Real*const u) : grid(g), Uinf(u)
  {
    MPI_Comm_rank(grid->getCartComm(),&rank);
  }

  ~IF3D_ObstacleFactory()
  {}

  std::vector<IF3D_ObstacleOperator * > create(ArgumentParser & parser);
};


#endif /* defined(__IF3D_ROCKS__IF3D_ObstacleFactory__) */
