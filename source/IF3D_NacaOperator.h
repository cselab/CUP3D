//
//  IF3D_CarlingFishOperator.h
//  IncompressibleFluids3D
//
//  Created by Wim van Rees on 4/15/13.
//
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
	IF3D_NacaOperator(FluidGridMPI * grid, ArgumentParser & parser);
	void _parseArguments(ArgumentParser & parser);
	void update(const int stepID, const Real t, const Real dt, const Real* Uinf) override;
	void computeVelocities(const Real Uinf[3]) override;
	void create(const int step_id,const Real time, const Real dt, const Real *Uinf) override;
   	void finalize(const int step_id,const Real time, const Real dt, const Real *Uinf) override;
};


#endif /* defined(__IncompressibleFluids3D__IF3D_CarlingFish__) */
