//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch) and Wim van Rees.
//

#ifndef __IncompressibleFluids3D__IF3D_StefanFishOperator__
#define __IncompressibleFluids3D__IF3D_StefanFishOperator__

#include "IF3D_FishOperator.h"

class IF3D_StefanFishOperator: public IF3D_FishOperator
{
protected:
public:
  IF3D_StefanFishOperator(FluidGridMPI*g, ArgumentParser&p, const Real*const u);
  void save(const int step_id, const double t, std::string filename = std::string()) override;
  void restart(const double t, std::string filename) override;

  #ifdef RL_LAYER
    void execute(const int iAgent, const double time, const vector<double>a) override;
  #endif
};

#endif /* defined(__IncompressibleFluids3D__IF3D_CarlingFish__) */
