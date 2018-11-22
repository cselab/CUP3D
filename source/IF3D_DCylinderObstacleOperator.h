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
            FluidGridMPI *g,
            ArgumentParser &p,
            const Real *u);

  void create(const int step_id,const double time, const double dt, const Real *Uinf) override;
};
