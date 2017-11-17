//
//  CubismUP_3D
//
//  Written by Guido Novati ( novatig@ethz.ch ).
//  This file started as an extension of code written by Wim van Rees
//  Copyright (c) 2017 ETHZ. All rights reserved.
//

#ifndef __IncompressibleFluids3D__IF3D_DeadFishOperator__
#define __IncompressibleFluids3D__IF3D_DeadFishOperator__

#include <cmath>
#include <array>
#include "IF3D_FishOperator.h"

class IF3D_DeadFishOperator: public IF3D_FishOperator
{
protected:
    double angVel_comp[3], transVel_comp[3], ext_pos[3]; // from computeVelocities
    double Ltow, Ttow, Atow, VelX, P0, Y0;
    bool bTilt;
    int ID;

public:
    IF3D_DeadFishOperator(FluidGridMPI*g, ArgumentParser&p, const Real*const u);
    void update(const int step_id, const double t, const double dt, const Real *Uinf) override;
    void computeVelocities(const Real* Uinf) override;
    void _parseArguments(ArgumentParser & parser);
  void save(const int step_id, const double t, std::string filename = std::string()) override;
  void restart(const double t, std::string filename = std::string()) override;
};


#endif /* defined(__IncompressibleFluids3D__IF3D_CarlingFish__) */
