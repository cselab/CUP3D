//
//  Cubism3D
//  Copyright (c) 2022 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//

#pragma once

#include "Fish.h"
#include "StefanFish.h"
#include "FishLibrary.h"
#include "FishShapes.h"
#include "ObstacleVector.h"

#include <Cubism/ArgumentParser.h>

#include <array>
#include <cmath>
#include <sstream>
#include <iomanip>

CubismUP_3D_NAMESPACE_BEGIN

class StefanFish: public Fish
{
public:
  StefanFish(SimulationData&s, cubism::ArgumentParser&p);

  //PID controller options
  bool bCorrectPosition ; //control yaw angle and position in x-y plane
  bool bCorrectPositionZ; //control pitch angle and position in z
  bool bCorrectRoll     ; //prevent rolling angle 
  Real origC[3]         ; //initial location for PID controller
  Real wyp              ; //weight for y-control
  Real wzp              ; //weight for z-control

  void create() override;
  virtual void computeVelocities() override;
  virtual void saveRestart( FILE * f ) override;
  virtual void loadRestart( FILE * f ) override;

  //Reinforcement Learning functions
  void act(const Real lTact, const std::vector<Real>& a) const; //get vector of actions that will be taken at time lTact
  std::vector<Real> state() const;                              //return vector of state
  Real getPhase(const Real time) const;
  Real getLearnTPeriod() const; //take actions once every 0.5*getLearnTPeriod() time
  //// Helpers for state function
  ssize_t holdingBlockID(const std::array<Real,3> pos) const;
  std::array<Real, 3> getShear(const std::array<Real,3> pSurf) const;
};


class CurvatureDefinedFishData : public FishMidlineData
{
 public:
  // stored past action for RL state:
  Real lastTact = 0;
  Real lastCurv = 0;
  Real oldrCurv = 0;
  // quantities needed to correctly control the speed of the midline maneuvers:
  Real periodPIDval = Tperiod;
  Real periodPIDdif = 0;
  bool TperiodPID = false;
  Real lastTime = 0;
  // quantities needed for rl:
  Real time0 = 0;
  Real timeshift = 0;

  // next scheduler is used to ramp-up the curvature from 0 during first period:
  Schedulers::ParameterSchedulerVector<6>    curvatureScheduler;
  // next scheduler is used for midline-bending control points for RL:
  Schedulers::ParameterSchedulerLearnWave<7> rlBendingScheduler;

  // pitching can be performed either by controling the midline torsion OR 
  // by performing a pitching motion by wrapping the fish around a cylinder
  bool control_torsion{false};
  // 
  // I.Torsion control parameters: torsion is a natural cubic spline passing through
  //   six points. Actions modify the torsion values directly.
  Schedulers::ParameterSchedulerVector<3>    torsionScheduler;
  std::array<Real,3> torsionValues          = {0,0,0};
  std::array<Real,3> torsionValues_previous = {0,0,0};
  Real Ttorsion_start = 0.0;

  Real  alpha     = 1;
  Real dalpha     = 0;
  Real  beta      = 0;
  Real dbeta      = 0;
  Real  gamma     = 0;
  Real dgamma     = 0;

  // next scheduler is used to ramp-up the period
  Schedulers::ParameterSchedulerScalar periodScheduler;
  Real current_period    = Tperiod;
  Real next_period       = Tperiod;
  Real transition_start  = 0.0;
  Real transition_duration = 0.1*Tperiod;

 protected:
  Real * const rK; //curvature kappa(s,t) of midline
  Real * const vK; //time derivative of curvature
  Real * const rC; //control parameters of curvature/bending
  Real * const vC; //control parameters of curvature/bending (=drC/dt)
  Real * const rB; //control parameters of curvature/bending
  Real * const vB; //control parameters of curvature/bending (=drB/dt)

  Real * const rT; //torsion tau(s,t) of midline
  Real * const vT; //time derivative of torsion
  Real * const rC_T; //control parameters of torsion
  Real * const vC_T; //control parameters of torsion (=drC_T/dt)
  Real * const rB_T; //control parameters of torsion
  Real * const vB_T; //control parameters of torsion (=drB_T/dt)

 public:
  CurvatureDefinedFishData(Real L, Real T, Real phi, Real _h, const Real _ampFac)
  : FishMidlineData(L, T, phi, _h, _ampFac),
    rK(_alloc(Nm)),vK(_alloc(Nm)), rC(_alloc(Nm)),vC(_alloc(Nm)), rB(_alloc(Nm)),vB(_alloc(Nm)),
    rT(_alloc(Nm)),vT(_alloc(Nm)), rC_T(_alloc(Nm)),vC_T(_alloc(Nm)), rB_T(_alloc(Nm)),vB_T(_alloc(Nm))
    {}

  void correctTailPeriod(const Real periodFac, const Real periodVel, const Real t, const Real dt)
  {
    assert(periodFac>0 && periodFac<2); // would be crazy

    const Real lastArg = (lastTime-time0)/periodPIDval + timeshift;
    time0 = lastTime;
    timeshift = lastArg;
    // so that new arg is only constant (prev arg) + dt / periodPIDval
    // with the new l_Tp:
    periodPIDval = Tperiod * periodFac;
    periodPIDdif = Tperiod * periodVel;
    lastTime = t;
    TperiodPID = true;
  }

  void execute(const Real time, const Real l_tnext, const std::vector<Real>& input) override;

  ~CurvatureDefinedFishData() override
  {
    _dealloc(rK); _dealloc(vK); _dealloc(rC);
    _dealloc(vC); _dealloc(rB); _dealloc(vB);
    _dealloc(rT); _dealloc(vT); _dealloc(rC_T);
    _dealloc(vC_T); _dealloc(rB_T); _dealloc(vB_T);
  }

  void computeMidline(const Real time, const Real dt) override;

  void performPitchingMotion(const Real time);

  void recomputeNormalVectors();

  //Actions available to StefanFish (for Reinforcement Learning)
  void action_curvature(const Real time, const Real l_tnext, const Real action)
  {
    rlBendingScheduler.Turn(action, l_tnext);
  }
  void action_period(const Real time, const Real l_tnext, const Real action)
  {
    if (TperiodPID) std::cout << "Warning: PID controller should not be used with RL." << std::endl;
    current_period = periodPIDval;
    next_period = Tperiod * (1 + action);
    transition_start = l_tnext;
  }
  void action_torsion(const Real time, const Real l_tnext, const Real * action)
  {
    for (int i = 0 ; i < 3 ; i++)
    {
      torsionValues_previous [i] = torsionValues[i];
      torsionValues[i] = action[i];
    }
    Ttorsion_start = time;
  }
  void action_torsion_pitching_radius(const Real time, const Real l_tnext, const Real action)
  {
    //Set torsion values so that a pitching motion is performed. This is similar to pitching for
    //a given cylinder radius R (as in action_pitching), but what is specified here is the 
    //midline torsion. The 'action' parameter is proportional to 1/R.
    const Real sq = 1.0/pow(2.0,0.5);
    const Real ar [3] = { action * ( norX[0     ] * 0.000 + norZ[0     ] * 1.000) ,
                          action * ( norX[Nm/2-1] * sq    + norZ[Nm/2-1] * sq   ) ,
                          action * ( norX[Nm-1  ] * 1.000 + norZ[Nm-1  ] * 0.000) };
    action_torsion(time, l_tnext, ar);
  }
};


CubismUP_3D_NAMESPACE_END
