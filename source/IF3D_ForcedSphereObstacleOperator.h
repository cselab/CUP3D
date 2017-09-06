//
//  IF3D_ForcedSphereObstacleOperator.h
//  IncompressibleFluids3D
//
//  Created by Wim van Rees on 8/20/13.
//
//

#ifndef __IncompressibleFluids3D__IF3D_ForcedSphereObstacleOperator__
#define __IncompressibleFluids3D__IF3D_ForcedSphereObstacleOperator__

#include "IF3D_SphereObstacleOperator.h"

#include <cmath>

class IF3D_ForcedSphereObstacleOperator: public IF3D_SphereObstacleOperator
{
  bool accel_decel = false;
  Real umax = 0, tmax = 1;
public:

	IF3D_ForcedSphereObstacleOperator(FluidGridMPI*grid, ArgumentParser&parser):
	IF3D_SphereObstacleOperator(grid,parser) {}

  // no need to compute velocities, are fixed
  void computeVelocities(const Real Uinf[3]) override
	{
		computeVelocities_forced(Uinf);
	}
  void _parseArguments(ArgumentParser & parser) override;
};

#endif /* defined(__IncompressibleFluids3D__IF3D_ForcedSphereObstacleOperator__) */
