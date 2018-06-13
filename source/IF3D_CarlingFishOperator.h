//
//  CubismUP_3D
//
//  Written by Guido Novati ( novatig@ethz.ch ).
//  This file started as an extension of code written by Wim van Rees
//  Copyright (c) 2017 ETHZ. All rights reserved.
//

#ifndef __IncompressibleFluids3D__IF3D_CarlingFishOperator__
#define __IncompressibleFluids3D__IF3D_CarlingFishOperator__

#include <cmath>
#include <array>
#include "IF3D_FishOperator.h"

class IF3D_CarlingFishOperator: public IF3D_FishOperator
{
public:
  IF3D_CarlingFishOperator(FluidGridMPI*g, ArgumentParser&p, const Real*const u);
  void execute(const int i,const double t,const vector<double>a) override;
  void computeForces(const int stepID, const double time, const double dt,
    const Real* Uinf, const double NU, const bool bDump) override;
};


#endif /* defined(__IncompressibleFluids3D__IF3D_CarlingFish__) */
