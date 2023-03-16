//
//  Cubism3D
//  Copyright (c) 2023 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//

#include "ExternalForcing.h"

CubismUP_3D_NAMESPACE_BEGIN
using namespace cubism;

void ExternalForcing::operator()(const double dt)
{
  sim.startProfiler("Forcing Kernel");
  const int dir = sim.BCy_flag==wall ? 1 : 2;
  const Real H = sim.extents[dir];
  const Real gradPdt = 8*sim.uMax_forced*sim.nu/H/H * dt;
  const int DIRECTION = 0;
  const std::vector<BlockInfo> &  velInfo = sim.velInfo();
  #pragma omp parallel for 
  for(size_t i=0; i<velInfo.size(); i++)
  {
    VectorBlock& v = *(VectorBlock*)velInfo[i].ptrBlock;
    for(int z=0; z<VectorBlock::sizeZ; ++z)
    for(int y=0; y<VectorBlock::sizeY; ++y)
    for(int x=0; x<VectorBlock::sizeX; ++x)
    {
        v(x,y,z).u[DIRECTION] += gradPdt;
    }
  }

  sim.stopProfiler();
}

CubismUP_3D_NAMESPACE_END
