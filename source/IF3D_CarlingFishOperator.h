//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#ifndef __IncompressibleFluids3D__IF3D_CarlingFishOperator__
#define __IncompressibleFluids3D__IF3D_CarlingFishOperator__

#include <cmath>
#include <array>
#include "IF3D_FishOperator.h"

class CarlingFishMidlineData;

class IF3D_CarlingFishOperator: public IF3D_FishOperator
{
  CarlingFishMidlineData* readHingeParams(ArgumentParser&p);
  CarlingFishMidlineData* readBurstCoastParams(ArgumentParser&p);
public:
  IF3D_CarlingFishOperator(FluidGridMPI*g, ArgumentParser&p, const Real*const u);
  void computeForces(const int stepID, const double time, const double dt,
    const Real* Uinf, const double NU, const bool bDump) override;

  #ifdef RL_LAYER
    void execute(const int i,const double t,const vector<double>a) override;
  #endif
};


#endif /* defined(__IncompressibleFluids3D__IF3D_CarlingFish__) */
