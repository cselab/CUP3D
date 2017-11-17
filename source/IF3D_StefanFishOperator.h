//
//  CubismUP_3D
//
//  Written by Guido Novati ( novatig@ethz.ch ).
//  This file started as an extension of code written by Wim van Rees
//  Copyright (c) 2017 ETHZ. All rights reserved.
//

#ifndef __IncompressibleFluids3D__IF3D_StefanFishOperator__
#define __IncompressibleFluids3D__IF3D_StefanFishOperator__

#include <cmath>
#include <array>
#include "IF3D_FishOperator.h"


class IF3D_StefanFishOperator: public IF3D_FishOperator
{
protected:
public:
  IF3D_StefanFishOperator(FluidGridMPI*g, ArgumentParser&p, const Real*const u);
  void _parseArguments(ArgumentParser & parser);
  void execute(const int iAgent, const double time, const vector<double>a) override;
	void save(const int step_id, const double t, std::string filename = std::string()) override;
  void restart(const double t, string filename) override;
};

#endif /* defined(__IncompressibleFluids3D__IF3D_CarlingFish__) */
