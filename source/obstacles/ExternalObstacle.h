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
#include "CubismUP_3D.h"
#include "obstacles/Obstacle.h"

CubismUP_3D_NAMESPACE_BEGIN

class ExternalObstacle : public Obstacle
{
public:
    cubismup3d::ExternalObstacleSettings settings;

    ExternalObstacle(
            SimulationData&s,
            const ObstacleArguments &args);
    ExternalObstacle(
            SimulationData&s, ArgumentParser &p)
        : ExternalObstacle(s, ObstacleArguments(s, p)) {}

    void computeVelocities() override;
    void create() override;
    void finalize() override;
};

CubismUP_3D_NAMESPACE_END
#endif // CubismUP_3D_ExternalObstacle_h
