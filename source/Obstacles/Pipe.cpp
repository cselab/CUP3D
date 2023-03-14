//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//

#include "Pipe.h"
#include "extra/ObstacleLibrary.h"

#include <Cubism/ArgumentParser.h>

CubismUP_3D_NAMESPACE_BEGIN
using namespace cubism;

namespace PipeObstacle
{
struct FillBlocks : FillBlocksBase<FillBlocks>
{
  const Real radius, halflength, h, safety = (2+SURFDH)*h;
  const Real position[3];
  const Real box[3][2] = {
    {(Real)position[0] - std::sqrt(2) / 2 * radius + safety, (Real)position[0] + std::sqrt(2) / 2 * radius - safety},
    {(Real)position[1] - std::sqrt(2) / 2 * radius + safety, (Real)position[1] + std::sqrt(2) / 2 * radius - safety},
    {(Real)position[2] - halflength-safety, (Real)position[2]+halflength+safety}
  };

  FillBlocks(const Real r, const Real halfl, const Real _h, const Real p[3]):
  radius(r), halflength(halfl), h(_h), position{p[0],p[1],p[2]} {}

  inline bool isTouching(const BlockInfo & info, const ScalarBlock&b) const
  {
    Real MINP[3], MAXP[3];
    info.pos(MINP, 0, 0, 0);
    info.pos(MAXP, ScalarBlock::sizeX-1, ScalarBlock::sizeY-1, ScalarBlock::sizeZ-1);
    const Real intersect[3][2] = {
        {std::max(MINP[0], box[0][0]), std::min(MAXP[0], box[0][1])},
        {std::max(MINP[1], box[1][0]), std::min(MAXP[1], box[1][1])},
        {std::max(MINP[2], box[2][0]), std::min(MAXP[2], box[2][1])}
    };
    return not ( intersect[0][1]-intersect[0][0]>0 &&
                 intersect[1][1]-intersect[1][0]>0 &&
                 intersect[2][1]-intersect[2][0]>0 );
  }

  inline Real signedDistance(const Real xo, const Real yo, const Real zo) const
  {
    const Real x = xo - position[0], y = yo - position[1], z = zo - position[2];
    const Real planeDist = radius - std::sqrt(x*x+y*y);
    const Real vertiDist = halflength - std::fabs(z);
    return -std::min(planeDist, vertiDist);
  }
};
}

Pipe::Pipe(
    SimulationData&s, ArgumentParser &p)
    : Obstacle(s, p), radius(.5 * length),
      halflength(p("-halflength").asDouble(.5 * sim.extents[2]))
{
  section = p("-section").asString("circular");
  accel = p("-accel").asBool(false);
  if(accel) {
    if(not bForcedInSimFrame[0]) {
      printf("Warning: Pipe was not set to be forced in x-dir, yet the accel pattern is active.\n");
    }
    umax = - p("-xvel").asDouble(0.0);
    vmax = - p("-yvel").asDouble(0.0);
    wmax = - p("-zvel").asDouble(0.0);
    tmax = p("-T").asDouble(1.0);
  }
  _init();
}

void Pipe::_init(void)
{
  if (sim.verbose) printf("Created Pipe with radius %f and halflength %f\n", radius, halflength);

  // D-cyl can float around the domain, but does not support rotation. TODO
  bBlockRotation[0] = true;
  bBlockRotation[1] = true;
  bBlockRotation[2] = true;
}


void Pipe::create()
{
  const Real h = sim.hmin;
  const PipeObstacle::FillBlocks kernel(radius, halflength, h, position);
  create_base<PipeObstacle::FillBlocks>(kernel);
}


void Pipe::computeVelocities()
{
  if(accel) {
    if(sim.time<tmax) transVel_imposed[0] = umax*sim.time/tmax;
    else
    {
       transVel_imposed[0] = umax;
       transVel_imposed[1] = vmax;
       transVel_imposed[2] = wmax;
    }
  }

  Obstacle::computeVelocities();
}

void Pipe::finalize()
{
  // this method allows any computation that requires the char function
  // to be computed. E.g. compute the effective center of mass or removing
  // momenta from udef
}

CubismUP_3D_NAMESPACE_END
