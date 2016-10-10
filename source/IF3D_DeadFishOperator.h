//
//  IF3D_CarlingFishOperator.h
//  IncompressibleFluids3D
//
//  Created by Wim van Rees on 4/15/13.
//
//

#ifndef __IncompressibleFluids3D__IF3D_CarlingFishOperator__
#define __IncompressibleFluids3D__IF3D_CarlingFishOperator__

#include <cmath>
#include <array>
#include "IF3D_FishOperator.h"

class IF3D_DeadFishOperator: public IF3D_FishOperator
{
protected:
    double angVel_comp[3], transVel_comp[3], ext_pos[3]; // from computeVelocities
    double Ltow, Ttow, Atow, VelX, P0;
    bool bTilt;
    int ID;

public:
    IF3D_DeadFishOperator(FluidGridMPI * grid, ArgumentParser & parser);
    void update(const int step_id, const Real t, const Real dt, const Real *Uinf) override;
    void computeVelocities(const Real* Uinf) override;
    void _parseArguments(ArgumentParser & parser);
	void save(const int step_id, const Real t, std::string filename = std::string()) override;
	void restart(const Real t, std::string filename = std::string()) override;
};


#endif /* defined(__IncompressibleFluids3D__IF3D_CarlingFish__) */
