//
//  CubismUP_3D
//
//  Written by Guido Novati ( novatig@ethz.ch ).
//  This file started as an extension of code written by Wim van Rees
//  Copyright (c) 2017 ETHZ. All rights reserved.
//

#include "IF3D_CarlingFishOperator.h"
#include "IF3D_FishLibrary.h"
#include "GenericOperator.h"

IF3D_CarlingFishOperator::IF3D_CarlingFishOperator(FluidGridMPI*g, ArgumentParser&p, const Real*const u)
: IF3D_FishOperator(g, p, u)
{
  _parseArguments(p);

  sr.updateInstant(position[0], absPos[0], position[1], absPos[1],
                    _2Dangle, transVel[0], transVel[1], angVel[2]);
}

void IF3D_CarlingFishOperator::_parseArguments(ArgumentParser & parser)
{
  IF3D_FishOperator::_parseArguments(parser);

  const int Nextension = NEXTDX*NPPEXT;// up to 3dx on each side (to get proper interpolation up to 2dx)
  const double target_Nm = TGTPPB*length/vInfo[0].h_gridpoint;
  const double dx_extension = (1./NEXTDX)*vInfo[0].h_gridpoint;
  const int Nm = (Nextension+1)*(int)std::ceil(target_Nm/(Nextension+1)) + 1;
  const double amplitude = parser("-amplitude").asDouble(0.1212121212121212);
  const bool _bBurst = parser("-BurstCoast").asBool(false);
  const bool bHinge = parser("-HingedFin").asBool(false);
  const bool bOptimizeHinge = parser("-OptimizeHingedFin").asBool(false);
  printf("CarlingFish: N:%d, L:%f, T:%f, phi:%f, dx_ext:%f, amplitude:%f\n",
          Nm,length,Tperiod,phaseShift,dx_extension,amplitude);
  fflush(0);

  if (_bBurst) {
    parser.set_strict_mode();
    const Real _tStart = parser("-tStartBC").asDouble();
    parser.unset_strict_mode();
    const string fburstpar = "burst_coast_carling_params.txt";
    myFish = new CarlingFishMidlineData(Nm, length, Tperiod, phaseShift, dx_extension, fburstpar, _tStart, amplitude);

  }else if (bHinge){
    parser.set_strict_mode();
    double sHinge = length*parser("-sHinge").asDouble();
    if(not bOptimizeHinge){
      double aHinge = parser("-AhingeDeg").asDouble();
      double phiHinge = parser("-phiHingeDeg").asDouble();
      myFish = new CarlingFishMidlineData(Nm, length, Tperiod, phaseShift, dx_extension,sHinge,aHinge,phiHinge, 0.0);
    }else{
      const bool equalHeight = length*parser("-equalHeight").asBool();
      //myFish = new CarlingFishMidlineData(Nm, length, Tperiod, phaseShift, dx_extension,sHinge, amplitude*0.625); // 0.625 necessary to have dx>0 when lambda>=0.5
      myFish = new CarlingFishMidlineData(Nm, length, Tperiod, phaseShift, dx_extension,sHinge, 0.0, equalHeight);
    }
    parser.unset_strict_mode();
  } else {
    myFish = new CarlingFishMidlineData(Nm, length, Tperiod, phaseShift, dx_extension, amplitude);
  }
}

void IF3D_CarlingFishOperator::execute(const int i, const double t, const vector<double>a)
{
  sr.resetAverage();
  sr.t_next_comm=1e6;
}

void IF3D_CarlingFishOperator::computeForces(const int stepID, const double time, const double dt, const Real* Uinf, const double NU, const bool bDump)
{
  IF3D_ObstacleOperator::computeForces(stepID, time, dt, Uinf, NU, bDump);

  // The torque at hinge2
  this->torqueZsection = sum[23];
}
