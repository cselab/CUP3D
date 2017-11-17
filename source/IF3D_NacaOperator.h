//
//  CubismUP_3D
//
//  Written by Guido Novati ( novatig@ethz.ch ).
//  This file started as an extension of code written by Wim van Rees
//  Copyright (c) 2017 ETHZ. All rights reserved.
//

#ifndef __IncompressibleFluids3D__IF3D_NacaOperator__
#define __IncompressibleFluids3D__IF3D_NacaOperator__

#include <cmath>
#include <array>
#include "IF3D_FishOperator.h"

class IF3D_NacaOperator: public IF3D_FishOperator
{
  double Apitch, Fpitch, Ppitch, Mpitch, Fheave, Aheave;
bool bCreated;
 public:
  IF3D_NacaOperator(FluidGridMPI*g, ArgumentParser&p, const Real*const u);
  void _parseArguments(ArgumentParser & parser);
  void update(const int stepID, const double t, const double dt, const Real* Uinf) override;
  void computeVelocities(const Real Uinf[3]) override;
  void create(const int step_id,const double time, const double dt, const Real *Uinf) override;
  void finalize(const int step_id,const double time, const double dt, const Real *Uinf) override;
};


#endif /* defined(__IncompressibleFluids3D__IF3D_CarlingFish__) */
