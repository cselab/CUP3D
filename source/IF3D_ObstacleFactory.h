//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch) and Wim van Rees.
//

#ifndef __IF3D_ROCKS__IF3D_ObstacleFactory__
#define __IF3D_ROCKS__IF3D_ObstacleFactory__

#include "IF3D_ObstacleOperator.h"

namespace cubism { struct ArgumentParser; }

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

  std::vector<IF3D_ObstacleOperator *> create(ArgumentParser &parser);
};


#endif /* defined(__IF3D_ROCKS__IF3D_ObstacleFactory__) */
