//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#ifndef CubismUP_3D_Sphere_h
#define CubismUP_3D_Sphere_h

#include "Obstacle.h"

CubismUP_3D_NAMESPACE_BEGIN


struct SphereArguments
{
  const Real radius;
  Real umax = 0;
  Real tmax = 1;
  //special case: startup with unif accel to umax in tmax, and then decel to 0
  bool accel_decel = false;
  bool bHemi = false;

  SphereArguments(Real R) : radius(R) {}
  SphereArguments(Real R, Real _umax, Real _tmax, bool _ad, bool _bHemi)
      : radius(R), umax(_umax), tmax(_tmax), accel_decel(_ad), bHemi(_bHemi) { }
};

class Sphere : public Obstacle, private SphereArguments
{
public:
  Sphere(SimulationData&s,cubism::ArgumentParser&p);
  void create() override;
  void finalize() override;
  void computeVelocities() override;
};

CubismUP_3D_NAMESPACE_END
#endif // CubismUP_3D_Sphere_h
