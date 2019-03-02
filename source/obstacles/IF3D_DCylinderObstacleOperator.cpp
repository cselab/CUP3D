//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#include "obstacles/IF3D_DCylinderObstacleOperator.h"
#include "obstacles/extra/IF3D_ObstacleLibrary.h"
#include "Cubism/ArgumentParser.h"

namespace DCylinderObstacle
{
struct FillBlocks : FillBlocksBase<FillBlocks>
{
  const Real radius, halflength, safety;
  const double position[3];
  const Real box[3][2] = {
    {position[0] - radius     - safety, position[0]              + safety},
    {position[1] - radius     - safety, position[1] + radius     + safety},
    {position[2] - halflength - safety, position[2] + halflength + safety}
  };

  FillBlocks(const Real r, const Real halfl, const Real h, const double p[3]):
  radius(r), halflength(halfl), safety(2*h), position{p[0],p[1],p[2]} {}

  bool _is_touching(const Real min_pos[3], const Real max_pos[3]) const
  {
    Real intersection[3][2] = {
        std::max(min_pos[0], box[0][0]), std::min(max_pos[0], box[0][1]),
        std::max(min_pos[1], box[1][0]), std::min(max_pos[1], box[1][1]),
        std::max(min_pos[2], box[2][0]), std::min(max_pos[2], box[2][1])
    };

    return
        intersection[0][1]-intersection[0][0]>0 &&
        intersection[1][1]-intersection[1][0]>0 &&
        intersection[2][1]-intersection[2][0]>0;
  }

  bool isTouching(const BlockInfo& info, const int buffer_dx = 0) const
  {
    Real min_pos[3], max_pos[3];
    info.pos(min_pos, 0,0,0);
    info.pos(max_pos, FluidBlock::sizeX-1, FluidBlock::sizeY-1, FluidBlock::sizeZ-1);
    for(int i=0;i<3;++i) {
      min_pos[i]-=buffer_dx*info.h_gridpoint;
      max_pos[i]+=buffer_dx*info.h_gridpoint;
    }
    return _is_touching(min_pos,max_pos);
  }

  inline Real signedDistance(const Real xo, const Real yo, const Real zo) const
  {
    const Real x = xo - position[0], y = yo - position[1], z = zo - position[2];
    const Real planeDist = std::min( -x, radius-std::sqrt(x*x+y*y) );
    const Real vertiDist = halflength - std::fabs(z);
    return std::min(planeDist, vertiDist);
  }
};
}

IF3D_DCylinderObstacleOperator::IF3D_DCylinderObstacleOperator(
    SimulationData&s,
    ArgumentParser &p)
    : IF3D_ObstacleOperator(s, p),
      radius(.5 * length),
      halflength(p("-halflength").asDouble(.5 * sim.extent[2]))
{
  _init();
}


IF3D_DCylinderObstacleOperator::IF3D_DCylinderObstacleOperator(
    SimulationData& s,
    ObstacleArguments &args,
    const double radius_,
    const double halflength_)
    : IF3D_ObstacleOperator(s, args), radius(radius_), halflength(halflength_)
{
  _init();
}

void IF3D_DCylinderObstacleOperator::_init(void) {
  printf("Created IF3D_DCylinderObstacleOperator with radius %f and halflength %f\n", radius, halflength);

  // D-cyl can float around the domain, but does not support rotation. TODO
  bBlockRotation[0] = true;
  bBlockRotation[1] = true;
  bBlockRotation[2] = true;
}


void IF3D_DCylinderObstacleOperator::create()
{
  const DCylinderObstacle::FillBlocks kernel(radius, halflength, vInfo[0].h_gridpoint, position);

  create_base<DCylinderObstacle::FillBlocks>(kernel);
}

void IF3D_DCylinderObstacleOperator::finalize()
{
  // this method allows any computation that requires the char function
  // to be computed. E.g. compute the effective center of mass or removing
  // momenta from udef
}
