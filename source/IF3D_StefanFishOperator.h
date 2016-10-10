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

namespace Fish
{
struct CurvatureDefinedFishData : FishMidlineData
{
protected:
	Schedulers::ParameterSchedulerVector<6> curvScheduler;
	Schedulers::ParameterSchedulerLearnWave<7> baseScheduler;
	Schedulers::ParameterSchedulerVector<6> adjustScheduler;
	Real * const rK;
	Real * const vK;
	Real * const rC;
	Real * const vC;
	Real * const rB;
	Real * const vB;
	Real * const rA;
	Real * const vA;
	Real l_Tp, time0, timeshift;

public:

	CurvatureDefinedFishData(const int Nm, const Real length, const Real Tperiod, const Real phaseShift, const Real dx_ext);
	void _correctTrajectory(const Real dtheta, const Real time, const Real dt) override
	void execute(const Real time, const Real l_tnext, const vector<Real>& input) override;
	~CurvatureDefinedFishData();
	void computeMidline(const Real time);
};
}

class IF3D_StefanFishOperator: public IF3D_FishOperator
{
protected:
	//bool randomActions, bSpiral, useLoadedActions;
    //vector<vector<Real>> loadedActions;
public:
    IF3D_StefanFishOperator(FluidGridMPI * grid, ArgumentParser & parser);
    void _parseArguments(ArgumentParser & parser);
    void execute(Communicator * comm, const int iAgent, const Real time) override;
};


#endif /* defined(__IncompressibleFluids3D__IF3D_CarlingFish__) */
