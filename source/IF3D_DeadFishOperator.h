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

class IF3D_CarlingFishOperator: public IF3D_FishOperator
{
public:
	
    IF3D_CarlingFishOperator(FluidGridMPI * grid, ArgumentParser & parser);
    ~IF3D_CarlingFishOperator();
	void save(const int step_id, const Real t, std::string filename = std::string()) override;
	void restart(const Real t, std::string filename = std::string()) override;
    void update(const int step_id, const Real t, const Real dt, const Real *Uinf) override;
    void getCenterOfMass(Real CM[3]) const override;
    void create(const int step_id,const Real time, const Real dt, const Real *Uinf) override;
    void _parseArguments(ArgumentParser & parser);
};


#endif /* defined(__IncompressibleFluids3D__IF3D_CarlingFish__) */
