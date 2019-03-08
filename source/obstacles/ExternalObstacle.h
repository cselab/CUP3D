//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Ivica Kicic (kicici@ethz.ch) in May 2018.
//

#ifndef CubismUP_3D_ExternalObstacle_h
#define CubismUP_3D_ExternalObstacle_h

/*
 * This obstacle can be used to insert obstacles whose shape and velocity is
 * defined by an external code. Intended to be used when CubismUP_3D used as a
 * library.
 */
#include "obstacles/Obstacle.h"
#include <functional>

CubismUP_3D_NAMESPACE_BEGIN

/*
 * Callbacks and other information for `IF3D_ExternalObstacleOperator`.
 *
 * This structure enables the user to define a custom obstacle.
 */
struct ExternalObstacleSettings
{
  typedef std::array<Real, 3> Point;
  typedef std::array<Real, 3> Velocity;

  /*
   * Check if given box is touching (intersecting) the object.
   *
   * False positives are allowed.
   */
  std::function<bool(Point low, Point high)> is_touching_fn;

  /*
   * Returns the signed distance to the object boundary.
   *
   * Positive values are to be returned for points inside the object,
   * negative for points outside of the object. Must be precise only close to
   * the obstacle surface.
   */
  std::function<Real(Point)> signed_distance_fn;

  /* Returns the local object velocity at the given location. */
  std::function<Velocity(Point)> velocity_fn;

  /* Returns the center-of-mass velocity of the object. */
  std::function<Point()> com_velocity_fn;

  /* Approx. length of the object. */
  Real length = 0.0;

  /* Object center location. (center of mass? probably not so important) */
  std::array<Real, 3> position{{(Real)0.5, (Real)0.5, (Real)0.5}};
};


class ExternalObstacle : public Obstacle
{
public:
  ExternalObstacleSettings settings;

  ExternalObstacle(SimulationData&s, const ObstacleArguments &args);
  ExternalObstacle(SimulationData&s, cubism::ArgumentParser &p)
      : ExternalObstacle(s, ObstacleArguments(s, p)) {}

  void computeVelocities() override;
  void create() override;
  void finalize() override;
};

CubismUP_3D_NAMESPACE_END
#endif // CubismUP_3D_ExternalObstacle_h
