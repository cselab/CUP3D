//
//  IF3D_CarlingFishOperator.cpp
//  IncompressibleFluids3D
//
//  Created by Wim van Rees on 4/15/13.
//
//

#include "IF3D_CarlingFishOperator.h"
#include "IF3D_FishLibrary.h"
#include "GenericOperator.h"

IF3D_CarlingFishOperator::IF3D_CarlingFishOperator(FluidGridMPI * grid, ArgumentParser & parser)
: IF3D_FishOperator(grid, parser)
{
	_parseArguments(parser);
	const Real target_Nm = TGTPPB*length/vInfo[0].h_gridpoint;
	const Real dx_extension = 0.25*vInfo[0].h_gridpoint;
	const int Nm = NPPSEG*(int)std::ceil(target_Nm/NPPSEG)+1;
	printf("%d %f %f %f %f\n",Nm,length,Tperiod,phaseShift,dx_extension);
	fflush(0);
	// multiple of NPPSEG: TODO why?
	myFish = new CarlingFishMidlineData(Nm, length, Tperiod, phaseShift, dx_extension);
}

IF3D_CarlingFishOperator::~IF3D_CarlingFishOperator()
{
	if(myFish not_eq nullptr) delete myFish;
}

void IF3D_CarlingFishOperator::_parseArguments(ArgumentParser & parser)
{
	IF3D_FishOperator::_parseArguments(parser);
}
