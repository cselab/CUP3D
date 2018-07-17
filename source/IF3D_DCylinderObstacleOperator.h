//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#pragma once

#include "IF3D_ObstacleOperator.h"

class IF3D_DCylinderObstacleOperator : public IF3D_ObstacleOperator
{
    const double radius;
    const double halflength;

public:
    IF3D_DCylinderObstacleOperator(
            FluidGridMPI * const g,
            ArgumentParser &p,
            const Real * const u) //const Real radius, const double position[3], const Real smoothing_length=-1):
        : IF3D_ObstacleOperator(g, p, u),
          radius(.5 * length),
          halflength(p("-halflength").asDouble(.5 * ext_Z))
    {
        printf("Created IF3D_DCylinderObstacleOperator with radius %f and halflength %f\n", radius, halflength);
        // D-cyl can float around the domain, but does not support rotation. TODO
        bBlockRotation[0] = true;
        bBlockRotation[1] = true;
        bBlockRotation[2] = true;
    }

  void create(const int step_id,const double time, const double dt, const Real *Uinf) override;
};
