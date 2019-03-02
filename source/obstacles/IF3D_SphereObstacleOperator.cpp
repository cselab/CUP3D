//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#include "obstacles/IF3D_SphereObstacleOperator.h"
#include "obstacles/extra/IF3D_ObstacleLibrary.h"
#include "Cubism/ArgumentParser.h"

namespace SphereObstacle
{
struct FillBlocks : FillBlocksBase<FillBlocks>
{
  const Real radius,safe_radius;
  const double sphere_position[3];
  const Real sphere_box[3][2] = {
    {sphere_position[0] - safe_radius, sphere_position[0] + safe_radius},
    {sphere_position[1] - safe_radius, sphere_position[1] + safe_radius},
    {sphere_position[2] - safe_radius, sphere_position[2] + safe_radius}
  };

  FillBlocks(const Real _radius, const Real max_dx, const double pos[3]):
    radius(_radius), safe_radius(radius+4*max_dx),
    sphere_position{pos[0],pos[1],pos[2]} { }

  inline bool is_touching(const Real min_pos[3], const Real max_pos[3]) const
  {
    Real intersection[3][2] = {
    std::max(min_pos[0],sphere_box[0][0]),std::min(max_pos[0],sphere_box[0][1]),
    std::max(min_pos[1],sphere_box[1][0]),std::min(max_pos[1],sphere_box[1][1]),
    std::max(min_pos[2],sphere_box[2][0]),std::min(max_pos[2],sphere_box[2][1])
    };
    return intersection[0][1]-intersection[0][0]>0 && intersection[1][1]-intersection[1][0]>0 && intersection[2][1]-intersection[2][0]>0;
  }

  inline bool isTouching(const BlockInfo& info, const int buffer_dx = 0) const
  {
    Real min_pos[3], max_pos[3];

    info.pos(min_pos, 0,0,0);
    info.pos(max_pos, FluidBlock::sizeX-1, FluidBlock::sizeY-1, FluidBlock::sizeZ-1);
    for(int i=0;i<3;++i) {
      min_pos[i]-=buffer_dx*info.h_gridpoint;
      max_pos[i]+=buffer_dx*info.h_gridpoint;
    }
    return is_touching(min_pos,max_pos);
  }

  inline Real signedDistance(const Real x, const Real y, const Real z) const
  {
    const Real dx = x-sphere_position[0];
    const Real dy = y-sphere_position[1];
    const Real dz = z-sphere_position[2];
    return radius - std::sqrt(dx*dx + dy*dy + dz*dz); // pos inside, neg outside
  }
};
}

IF3D_SphereObstacleOperator::IF3D_SphereObstacleOperator(
    SimulationData& s, ArgumentParser& p )
    : IF3D_ObstacleOperator(s, p), radius(0.5*length)
{
  accel_decel = p("-accel").asBool(false);
  if(accel_decel) {
    if(not bForcedInSimFrame[0]) {
      printf("Warning: sphere was not set to be forced in x-dir, yet the accel_decel pattern is active.\n");
    }
    umax = p("-xvel").asDouble(0.0);
    tmax = p("-T").asDouble(1.);
  }
}


IF3D_SphereObstacleOperator::IF3D_SphereObstacleOperator(
    SimulationData& s,
    ObstacleArguments &args,
    const double R)
    : IF3D_ObstacleOperator(s, args), radius(R) { }


IF3D_SphereObstacleOperator::IF3D_SphereObstacleOperator(
    SimulationData& s,
    ObstacleArguments &args,
    const double R,
    const double _umax,
    const double _tmax)
    : IF3D_ObstacleOperator(s, args), radius(R), umax(_umax), tmax(_tmax) { }


void IF3D_SphereObstacleOperator::create()
{
  const SphereObstacle::FillBlocks K(radius, vInfo[0].h_gridpoint, position);

  create_base<PlateObstacle::FillBlocks>(K);
}

void IF3D_SphereObstacleOperator::finalize()
{
  // this method allows any computation that requires the char function
  // to be computed. E.g. compute the effective center of mass or removing
  // momenta from udef
}


void IF3D_SphereObstacleOperator::computeVelocities()
{
  IF3D_ObstacleOperator::computeVelocities();

  if(accel_decel) {
    if(sim.time<tmax) transVel[0] = umax*sim.time/tmax;
    else if (sim.time<2*tmax) transVel[0] = umax*(2*tmax-sim.time)/tmax;
    else transVel[0] = 0;
  }
}
