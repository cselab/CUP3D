//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#include "obstacles/IF3D_CarlingFishOperator.h"
#include "obstacles/extra/IF3D_FishLibrary.h"

#include "Cubism/ArgumentParser.h"

class CarlingFishMidlineData : public FishMidlineData
{
 public:
  bool quadraticAmplitude = true;
 protected:
  const double carlingAmp;
  static constexpr double carlingInv = 0.03125;

  const double quadraticFactor; // Should be set to 0.1, which gives peak-to-peak amp of 0.2L (this is physically observed in most fish species)

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
    // Maertens et al. JFM 2017:
    return quadraticFactor*(length -.825*(s-length) +1.625*(s*s/length-length));
    //return s*s*quadraticFactor/length;
  }
  inline Real getLinAmp(const Real s) const
  {
    return carlingAmp * (s + carlingInv*length);
  }

  inline Real getArg(const Real s,const Real t) const
  {
    return 2.0*M_PI*(s/(waveLength*length) - t/Tperiod + phaseShift);
  }

  virtual Real midline(const Real s, const Real t) const;

  virtual Real midlineVel(const Real s, const Real t) const;

  // This needed only during burstCoast
  std::pair<double, double> cubicHermite(const double f1, const double f2, const double x){
    const double a =  2*(f1-f2);
    const double b = -3*(f1-f2);
    const double retVal = a*x*x*x + b*x*x + f1;
    const double deriv = 3*a*x*x + 2*b*x;
    return std::make_pair(retVal, deriv);
  }

  virtual void _computeMidlineCoordinates(const Real time);

  virtual void _computeMidlineVelocities(const Real time);

 public:
  // L=length, T=period, phi=phase shift, _h=grid size, A=amplitude modulation
  CarlingFishMidlineData(double L, double T, double phi, double _h, double A) :
  FishMidlineData(L,T,phi,_h,A),carlingAmp(.1212121212*A),quadraticFactor(.1*A)
  {
    // FinSize has now been updated with value read from text file. Recompute heights to over-write with updated values
    //printf("Overwriting default tail-fin size for Plain Carling:\n");
    //_computeWidthsHeights();
  }

  virtual void computeMidline(const double t, const double dt) override;
};

#include "obstacles/extra/IF3D_CarlingFishOperator_extra.h"

void CarlingFishMidlineData::computeMidline(const double t,const double dt)
{
  _computeMidlineCoordinates(t);
  _computeMidlineVelocities(t);
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

Real CarlingFishMidlineData::midline(const Real s, const Real t) const
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

Real CarlingFishMidlineData::midlineVel(const Real s, const Real t) const
{
  const Real arg = getArg(s, t);
  const Real dArg = -2*M_PI/Tperiod;

  double velCurrent;
  if(quadraticAmplitude) {
    velCurrent = getQuadAmp(s) * dArg * std::cos(arg);
  }else{
    velCurrent =  getLinAmp(s) * dArg * std::cos(arg);
  }

  return velCurrent;
}

void CarlingFishMidlineData::_computeMidlineCoordinates(const Real t)
{
  const Real rampFac = rampFactorSine(t, Tperiod);
  rX[0] = 0.0;
  rY[0] = rampFac*midline(rS[0], t);

  for(int i=1;i<Nm;++i)
  {
    rY[i] = rampFac*midline(rS[i], t);
    const Real dy = rY[i]-rY[i-1], ds = rS[i]-rS[i-1];
    Real dx = std::sqrt(ds*ds-dy*dy);
    assert(dx>0);
    rX[i] = rX[i-1] + dx;
  }
}

void CarlingFishMidlineData::_computeMidlineVelocities(const Real t)
{
  const Real rampFac =    rampFactorSine(t, Tperiod);
  const Real rampFacVel = rampFactorVelSine(t, Tperiod);

  vX[0] = 0.0; //rX[0] is constant
  vY[0] = rampFac*midlineVel(rS[0],t) + rampFacVel*midline(rS[0],t);

  for(int i=1; i<Nm; ++i)
  {
    vY[i]=rampFac*midlineVel(rS[i],t) + rampFacVel*midline(rS[i],t);
    const Real dy = rY[i]-rY[i-1], dx = rX[i]-rX[i-1], dVy = vY[i]-vY[i-1];
    assert(dx>0);
    vX[i] = vX[i-1] - dy/dx *dVy; // use ds^2 = dx^2+dy^2 -> ddx = -dy/dx*ddy
  }
}

IF3D_CarlingFishOperator::IF3D_CarlingFishOperator(SimulationData&s,
  ArgumentParser&p) : IF3D_FishOperator(s, p)
{
  // _ampFac=0.0 for towed fish :
  const double ampFac = p("-amplitudeFactor").asDouble(1.0);
  const bool bQuadratic = p("-bQuadratic").asBool(true);
  const bool bBurst = p("-BurstCoast").asBool(false);
  const bool bHinge = p("-HingedFin").asBool(false);
  if(bBurst && bHinge) {
    printf("Pick either hinge or burst and coast!\n"); fflush(0);
    MPI_Abort(MPI_COMM_WORLD,1);
  }
  if(bBurst || bHinge) printf("WARNING: UNTESTED!!!\n");

  CarlingFishMidlineData* localFish = nullptr; //could be class var if needed
  if(bBurst) localFish = readBurstCoastParams(p);
  else
  if(bHinge) localFish = readHingeParams(p);
  else
  localFish = new CarlingFishMidlineData(length, Tperiod, phaseShift,
    vInfo[0].h_gridpoint, ampFac);

  // generic copy for base class:
  assert( myFish == nullptr );
  myFish = (FishMidlineData*) localFish;

  localFish->quadraticAmplitude = bQuadratic;
  std::string heightName = p("-heightProfile").asString("baseline");
  std::string  widthName = p( "-widthProfile").asString("baseline");
  MidlineShapes::computeWidthsHeights(heightName, widthName, length,
    myFish->rS, myFish->height, myFish->width, myFish->Nm, sim.rank);

  if(!sim.rank)
    printf("CarlingFish: N:%d, L:%f, T:%f, phi:%f, amplitude:%f\n",
        myFish->Nm, length, Tperiod, phaseShift, ampFac);

  #ifdef RL_LAYER
    sr = StateReward(length, Tperiod);
    sr.parseArguments(p);
    sr.updateInstant(position[0], absPos[0], position[1], absPos[1],
                      _2Dangle, transVel[0], transVel[1], angVel[2]);
  #endif
}


#ifdef RL_LAYER

void IF3D_CarlingFishOperator::execute(const int i, const double t, const std::vector<double>a)
{
  sr.resetAverage();
  sr.t_next_comm=1e6;
}

#endif
