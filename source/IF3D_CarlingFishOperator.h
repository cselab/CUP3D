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
    void _parseArguments(ArgumentParser & parser);
    void execute(Communicator * comm, const int iAgent, const Real time, const int iLabel) override;
};


#endif /* defined(__IncompressibleFluids3D__IF3D_CarlingFish__) */
