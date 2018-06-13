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


class CarlingFishMidlineData : public FishMidlineData
{
 public:
  bool quadraticAmplitude = true;
 protected:
  const double carlingAmp;
  static constexpr double carlingInv = 0.03125;

  static constexpr double quadraticFactor = 0.1;

  Real aParabola = 0;
  Real bParabola = 0;
  Real cParabola = 0;

  inline Real rampFactorSine(const Real t, const Real T) const
  {
    return (t<T ? std::sin(0.5*M_PI*t/T) : 1.0);
  }

  inline Real rampFactorVelSine(const Real t, const Real T) const
  {
    return (t<T ? 0.5*M_PI/T * std::cos(0.5*M_PI*t/T) : 0.0);
  }

  inline Real getQuadAmp(const Real s) const
  {
    return s*s*quadraticFactor/length;
  }
  inline Real getLinAmp(const Real s) const
  {
    return carlingAmp * (s + carlingInv*length);
  }

  inline Real getArg(const Real s,const Real t) const
  {
    return 2.0*M_PI*(s/(waveLength*length) - t/Tperiod + phaseShift);
  }

  //inline void computeParabolaParams(const Real yLeft, const Real yPrimeLeft, const Real yPrimeRight) const
  //{
  //  aParabola = (yPrimeLeft-yPrimeRight)/(2*(sLeft-sRight));
  //  bParabola = yPrimeRight - 2*aParabola*sRight;
  //  cParabola = yLeft - aParabola*sLeft*sLeft - bParabola*sLeft;
  //}

  Real getJointParabola(const Real s, const Real L) const
  {
    return aParabola*s*s + bParabola*s + cParabola;
  }

  Real midline(const Real s, const Real t) const
  {
    const Real arg = getArg(s, t);

    double yCurrent;
    if(quadraticAmplitude){
      yCurrent = getQuadAmp(s) * std::sin(arg);
    } else {
      yCurrent =  getLinAmp(s) * std::sin(arg);
    }

    return yCurrent;
  }

  inline Real midlineVel(const Real s, const Real t) const
  {
    const Real arg = getArg(s, t);
    const Real dArg = -2*M_PI/Tperiod;

    double velCurrent;
    if(quadraticAmplitude) {
      velCurrent = getQuadAmp(s) * dArg * std::cos(arg);
    }else{
      velCurrent = -getLinAmp(s) * dArg * std::cos(arg);
    }

    return velCurrent;
  }

  std::pair<double, double> cubicHermite(const double f1, const double f2, const double x){
    const double a =  2*(f1-f2);
    const double b = -3*(f1-f2);
    const double retVal = a*x*x*x + b*x*x + f1;
    const double deriv = 3*a*x*x + 2*b*x;
    return std::make_pair(retVal, deriv);
  }

  void _computeMidlineCoordinates(const Real time)
  {
    const Real rampFac = rampFactorSine(time, Tperiod);
    rX[0] = 0.0;
    rY[0] = rampFac*midline(rS[0], time);

    for(int i=1;i<Nm;++i)
    {
      rY[i]=rampFac*midline(rS[i], time);
      const Real dy = rY[i]-rY[i-1];
      const Real ds = rS[i]-rS[i-1];
      Real dx = std::sqrt(ds*ds-dy*dy);
      assert(dx>0);
      rX[i] = rX[i-1] + dx;
    }
  }

  void _computeMidlineVelocities(const Real time)
  {
    const Real rampFac =    rampFactorSine(time, Tperiod);
    const Real rampFacVel = rampFactorVelSine(time, Tperiod);

    vX[0] = 0.0; //rX[0] is constant
    vY[0] = rampFac*midlineVel(rS[0],time) + rampFacVel*midline(rS[0],time);

    for(int i=1;i<Nm;++i)
    {
      vY[i]=rampFac*midlineVel(rS[i],time) + rampFacVel*midline(rS[i],time);
      const Real dy  = rY[i]-rY[i-1];
      const Real dx  = rX[i]-rX[i-1];
      const Real dVy = vY[i]-vY[i-1];
      assert(dx>0);
      vX[i] = vX[i-1] - dy/dx *dVy; // use ds^2 = dx^2+dy^2 -> ddx = -dy/dx*ddy
    }
  }

 public:
  CarlingFishMidlineData(double L, double T, double phi, double _h, double A)
    : FishMidlineData(L, T, phi, _h), carlingAmp(A)
  {
    // FinSize has now been updated with value read from text file. Recompute heights to over-write with updated values
    //printf("Overwriting default tail-fin size for Plain Carling:\n");
    //_computeWidthsHeights();
  }

  void computeMidline(const Real time, const Real dt) override
  {
    _computeMidlineCoordinates(time);
    _computeMidlineVelocities(time);
    _computeMidlineNormals();
    #if 0
      #warning USED MPI COMM WORLD
      // we dump the profile
      int rank;
      MPI_Comm_rank(MPI_COMM_WORLD,&rank);
      if (rank!=0) return;
      FILE * f = fopen("fish_profile","w");
      for(int i=0;i<Nm;++i) fprintf(f,"%d %g %g %g %g %g %g %g %g %g\n",
       i,rS[i],rX[i],rY[i],norX[i],norY[i],vX[i],vY[i], vNorX[i],vNorY[i]);
      fclose(f); printf("Dumped midline\n");
    #endif
  }
};

IF3D_CarlingFishOperator::IF3D_CarlingFishOperator(FluidGridMPI*g,
  ArgumentParser&p, const Real*const u) : IF3D_FishOperator(g, p, u)
{
  sr = StateReward(length, Tperiod);
  sr.parseArguments(p);

  const double amplitude = p("-amplitude").asDouble(0.1212121212121212);
  myFish = new CarlingFishMidlineData(length, Tperiod, phaseShift,
    vInfo[0].h_gridpoint, amplitude);
  string heightName = p("-heightProfile").asString("baseline");
  string  widthName = p( "-widthProfile").asString("baseline");
  MidlineShapes::computeWidthsHeights(heightName, widthName, length,
    myFish->rS, myFish->height, myFish->width, myFish->Nm, rank);

  if(!rank)
    printf("CarlingFish: N:%d, L:%f, T:%f, phi:%f, amplitude:%f\n",
        myFish->Nm,length,Tperiod,phaseShift,amplitude);
            printf("%d %f %f %f\n",myFish->Nm, length, Tperiod, phaseShift);

  sr.updateInstant(position[0], absPos[0], position[1], absPos[1],
                    _2Dangle, transVel[0], transVel[1], angVel[2]);
}

void IF3D_CarlingFishOperator::execute(const int i, const double t, const vector<double>a)
{
  sr.resetAverage();
  sr.t_next_comm=1e6;
}

void IF3D_CarlingFishOperator::computeForces(const int stepID, const double time, const double dt, const Real* Uinf, const double NU, const bool bDump)
{
  IF3D_ObstacleOperator::computeForces(stepID, time, dt, Uinf, NU, bDump);
}
