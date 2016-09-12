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

#include "IF3D_Headers.h"
#include "IF3D_Types.h"
#include "IF3D_DeformingObstacleOperator.h"

class IF3D_CarlingFishOperator: public IF3D_ObstacleOperator
{
	Fish::CarlingFishMidlineData * myFish;
    Real length, Tperiod, phaseShift, phase, sim_time, sim_dt;
    
    double volume_internal, J_internal, CoM_internal[2], vCoM_internal[2], theta_internal, angvel_internal, angvel_internal_prev, CoM_interpolated[3];
    bool randomStart, randomActions, bSpiral, useLoadedActions, bCorrectTrajectory;
    //mt19937 * gen; //TODO... what to do? shared seed?
    Real Tstartlearn, GoalDX, new_curv, old_curv, new_Tp;
    vector<vector<Real>> loadedActions;
    int  nActions;

    //Real signLastTurn;
    //bool bObstacleBlocksFilled;
    //int nTurnActions, nPauseActions;
    //const bool bFixToPlanar; // restrict translation to xy plane, and rotation to z-axis

public:
	
    IF3D_CarlingFishOperator(FluidGridMPI * grid, ArgumentParser & parser)
    : IF3D_ObstacleOperator(grid, parser), theta_internal(0.0), angvel_internal(0.0), sim_time(0.0), sim_dt(0.0)
	{
        volume=0;
        for(int i=0;i<3;i++) transVel[i]=0;
        for(int i=0;i<3;i++) angVel[i]=0;
        for(int i=0;i<6;i++) J[i]=0;
        _parseArguments(parser);
#ifndef NDEBUG
        const double q_length=std::sqrt(quaternion[0]*quaternion[0]+quaternion[1]*quaternion[1]+quaternion[2]*quaternion[2]+quaternion[3]*quaternion[3]);
        assert(std::abs(q_length-1.0) < 5*std::numeric_limits<Real>::epsilon());
        assert(smoothing_length <= 4.0*(1./B::sizeX)*pow(0.5,grid.getCurrentMaxLevel()));
#endif
        const Real target_Nm = 2.0*length/vInfo[0].h_gridpoint;
        // multiple of 100: TODO why?
        const int Nm = 100*(int)std::ceil(target_Nm/100) + 1;
        const Real dx_extension = 0.25*vInfo[0].h_gridpoint;
        const int Nextension = 12;// up to 3dx on each side (to get proper interpolation up to 2dx)
        myFish = new Fish::CarlingFishMidlineData(Nm, sim_time, length, Tperiod, phaseShift, std::make_pair(Nextension,dx_extension));
    }
    
    ~IF3D_CarlingFishOperator();
	void save(const double t, std::string filename = std::string()) override;
	void restart(const double t, std::string filename = std::string()) override;
    void update(const double t, const double dt) override;
    void getCenterOfMass(double CM[3]) const override;

    void _parseArguments(ArgumentParser & parser);
};


#endif /* defined(__IncompressibleFluids3D__IF3D_CarlingFish__) */
