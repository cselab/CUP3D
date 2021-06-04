#pragma once

#include "Obstacle.h"
#include <functional>
#include <utility>

CubismUP_3D_NAMESPACE_BEGIN

using MovingSphereFunc = std::function<std::array<Real, 6>(Real time)>;

/// Create a sphere whose center position and COM velocity are governed by a
/// given function. This obstacle cannot be generated via the obstacle factory.
/// The function should return an array of 6 elements (position and velocity).
std::shared_ptr<cubismup3d::Obstacle> createMovingSphere(
    SimulationData &s,
    const ObstacleArguments &args,
    MovingSphereFunc func,
    Real radius);

CubismUP_3D_NAMESPACE_END
