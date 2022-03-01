//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#include "CarlingFish.h"
#include "FishLibrary.h"
#include "FishShapes.h"

#include <Cubism/ArgumentParser.h>

CubismUP_3D_NAMESPACE_BEGIN
using namespace cubism;

class CarlingFishMidlineData : public FishMidlineData
{
 public:
  bool quadraticAmplitude = false;
 protected:
  const double carlingAmp;
  static constexpr double carlingInv = 0.03125;

  const double quadraticFactor; // Should be set to 0.1, which gives peak-to-peak amp of 0.2L (this is physically observed in most fish species)

  inline Real rampFactorSine(const Real t, const Real T) const
  {
    //return (t<T ? ( 1 - std::cos(M_PI*t/T) )/2 : 1.0);
    return (t<T ? std::sin(0.5*M_PI*t/T) : 1.0);
  }

  inline Real rampFactorVelSine(const Real t, const Real T) const
  {
    //return (t<T ? 0.5*M_PI/T * std::sin(M_PI*t/T) : 0.0);
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

  template<bool bQuadratic>
  void _computeMidlinePosVel(const Real t)
  {
    const Real rampFac = rampFactorSine(t, Tperiod), dArg = -2*M_PI/Tperiod;
    const Real rampFacVel = rampFactorVelSine(t, Tperiod);
    {
      const Real arg = getArg(rS[0], t);
      const Real cosa = std::cos(arg), sina = std::sin(arg);
      const Real amp = bQuadratic? getQuadAmp(rS[0]) : getLinAmp(rS[0]);
      const Real Y = sina * amp, VY = cosa * dArg * amp;
      rX[0] = 0.0; vX[0] = 0.0; //rX[0] is constant
      rY[0] = rampFac*Y;
      vY[0] = rampFac*VY + rampFacVel*Y;
      rZ[0] = 0.0; vZ[0] = 0.0;
    }
    for(int i=1; i<Nm; ++i)
    {
      const Real arg = getArg(rS[i], t);
      const Real cosa = std::cos(arg), sina = std::sin(arg);
      const Real amp = bQuadratic? getQuadAmp(rS[i]) : getLinAmp(rS[i]);
      const Real Y = sina * amp, VY = cosa * dArg * amp;
      rY[i] = rampFac*Y;
      vY[i] = rampFac*VY + rampFacVel*Y;
      const Real dy = rY[i]-rY[i-1], ds = rS[i]-rS[i-1], dVy = vY[i]-vY[i-1];
      const Real dx = std::sqrt(ds*ds-dy*dy);
      assert(dx>0);
      rX[i] = rX[i-1] + dx;
      vX[i] = vX[i-1] - dy/dx *dVy; // use ds^2 = dx^2+dy^2 -> ddx = -dy/dx*ddy
      rZ[i] = 0.0;
      vZ[i] = 0.0;
    }
  }
};

void CarlingFishMidlineData::computeMidline(const double t,const double dt)
{
  if(quadraticAmplitude) _computeMidlinePosVel<true >(t);
  else                   _computeMidlinePosVel<false>(t);

  #pragma omp parallel for schedule(static)
  for(int i=0; i<Nm-1; i++) {
    const double ds = rS[i+1]-rS[i];
    const double tX = rX[i+1]-rX[i];
    const double tY = rY[i+1]-rY[i];
    const double tVX = vX[i+1]-vX[i];
    const double tVY = vY[i+1]-vY[i];
    norX[i] = -tY/ds;
    norY[i] =  tX/ds;
    norZ[i] =  0.0;
    vNorX[i] = -tVY/ds;
    vNorY[i] =  tVX/ds;
    vNorZ[i] = 0.0;
    binX[i] =  0.0;
    binY[i] =  0.0;
    binZ[i] =  1.0;
    vBinX[i] = 0.0;
    vBinY[i] = 0.0;
    vBinZ[i] = 0.0;
  }
  norX[Nm-1] = norX[Nm-2];
  norY[Nm-1] = norY[Nm-2];
  norZ[Nm-1] = norZ[Nm-2];
  vNorX[Nm-1] = vNorX[Nm-2];
  vNorY[Nm-1] = vNorY[Nm-2];
  vNorZ[Nm-1] = vNorZ[Nm-2];
  binX[Nm-1] = binX[Nm-2];
  binY[Nm-1] = binY[Nm-2];
  binZ[Nm-1] = binZ[Nm-2];
  vBinX[Nm-1] = vBinX[Nm-2];
  vBinY[Nm-1] = vBinY[Nm-2];
  vBinZ[Nm-1] = vBinZ[Nm-2];
}

CarlingFish::CarlingFish(SimulationData&s, ArgumentParser&p) : Fish(s, p)
{
  // _ampFac=0.0 for towed fish :
  const double ampFac = p("-amplitudeFactor").asDouble(1.0);
  const bool bQuadratic = p("-bQuadratic").asBool(false);
  const double Tperiod = p("-T").asDouble(1.0);
  const double phaseShift = p("-phi").asDouble(0.0);

  CarlingFishMidlineData* localFish = new CarlingFishMidlineData(length, Tperiod, phaseShift,
    sim.hmin, ampFac);

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
}

CubismUP_3D_NAMESPACE_END
