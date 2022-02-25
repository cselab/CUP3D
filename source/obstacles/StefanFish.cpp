//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch) and Wim van Rees.
//

#include "StefanFish.h"
#include "FishLibrary.h"
#include "FishShapes.h"

#include <Cubism/ArgumentParser.h>

#include <array>
#include <cmath>
#include <sstream>

CubismUP_3D_NAMESPACE_BEGIN
using namespace cubism;

class CurvatureDefinedFishData : public FishMidlineData
{
 public:
  // PID controller of body curvature:
  Real curv_PID_fac = 0;
  Real curv_PID_dif = 0;
  // exponential averages:
  Real avgDeltaY = 0;
  Real avgDangle = 0;
  Real avgAngVel = 0;
  // stored past action for RL state:
  Real lastTact = 0;
  Real lastCurv = 0;
  Real oldrCurv = 0;
  // quantities needed to correctly control the speed of the midline maneuvers:
  Real periodPIDval = Tperiod;
  Real periodPIDdif = 0;
  bool TperiodPID = false;
  // quantities needed for rl:
  Real time0 = 0;
  Real timeshift = 0;
  // aux quantities for PID controllers:
  Real lastTime = 0;
  Real lastAvel = 0;

  // next scheduler is used to ramp-up the curvature from 0 during first period:
  Schedulers::ParameterSchedulerVector<6>    curvatureScheduler;
  // next scheduler is used for midline-bending control points for RL:
  Schedulers::ParameterSchedulerLearnWave<7> rlBendingScheduler;

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

  CurvatureDefinedFishData(double L, double T, double phi, double _h, const double _ampFac)
  : FishMidlineData(L, T, phi, _h, _ampFac),
    rK(_alloc(Nm)),vK(_alloc(Nm)), rC(_alloc(Nm)),vC(_alloc(Nm)), rB(_alloc(Nm)),vB(_alloc(Nm)),
    rT(_alloc(Nm)),vT(_alloc(Nm)), rC_T(_alloc(Nm)),vC_T(_alloc(Nm)), rB_T(_alloc(Nm)),vB_T(_alloc(Nm)) { }

  void correctTrajectory(const Real dtheta, const Real vtheta,
                         const Real t, const Real dt)
  {
    curv_PID_fac = dtheta;
    curv_PID_dif = vtheta;
  }

  void correctTailPeriod(const Real periodFac,const Real periodVel,
                         const Real t, const Real dt)
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

  void execute(const double time, const double l_tnext, const std::vector<double>& input) override;

  ~CurvatureDefinedFishData() override
  {
    _dealloc(rK); _dealloc(vK); _dealloc(rC);
    _dealloc(vC); _dealloc(rB); _dealloc(vB);
    _dealloc(rT); _dealloc(vT); _dealloc(rC_T);
    _dealloc(vC_T); _dealloc(rB_T); _dealloc(vB_T);
  }

  void computeMidline(const double time, const double dt) override;
};

void CurvatureDefinedFishData::execute(const double time, const double l_tnext,
                                       const std::vector<double>& input)
{
  if (input.size()>1) {
    // printf("Turning by %g at tsim %g with tact %g.\n", input[0], time, l_tnext);
    rlBendingScheduler.Turn(input[0], l_tnext);
    //first, shift time to  previous turn node
    timeshift += (l_tnext-time0)/periodPIDval;
    time0 = l_tnext;
    periodPIDval = Tperiod*(1.+input[1]);
    periodPIDdif = 0;
  } else if (input.size()>0) {
    // printf("Turning by %g at tsim %g with tact %g.\n", input[0], time, l_tnext);
    rlBendingScheduler.Turn(input[0], l_tnext);
  }
}

void CurvatureDefinedFishData::computeMidline(const double t, const double dt)
{
  const std::array<Real ,6> curvaturePoints = { (Real)0, (Real).15*length,
    (Real).4*length, (Real).65*length, (Real).9*length, length
  };
  const std::array<Real,7> bendPoints = {(Real)-.5, (Real)-.25, (Real)0,
    (Real).25, (Real).5, (Real).75, (Real)1};
  const std::array<Real ,6> curvatureValues = {
      (Real)0.82014/length, (Real)1.46515/length, (Real)2.57136/length,
      (Real)3.75425/length, (Real)5.09147/length, (Real)5.70449/length
  };

  #if 1 // ramp-up over Tperiod
  const std::array<Real,6> curvatureZeros = std::array<Real, 6>();
  curvatureScheduler.transition(0,0,Tperiod,curvatureZeros ,curvatureValues);
  #else // no rampup for debug
  curvatureScheduler.transition(t,0,Tperiod,curvatureValues,curvatureValues);
  #endif

  // query the schedulers for current values
  curvatureScheduler.gimmeValues(t,                curvaturePoints,Nm,rS,rC,vC);
  rlBendingScheduler.gimmeValues(t,periodPIDval,length, bendPoints,Nm,rS,rB,vB);

  // next term takes into account the derivative of periodPIDval in darg:
  const Real diffT = TperiodPID? 1 - (t-time0)*periodPIDdif/periodPIDval : 1;
  // time derivative of arg:
  const Real darg = 2*M_PI/periodPIDval * diffT;
  const Real arg0 = 2*M_PI*((t-time0)/periodPIDval +timeshift) +M_PI*phaseShift;

  #pragma omp parallel for schedule(static)
  for(int i=0; i<Nm; ++i) {
    const Real arg = arg0 - 2*M_PI*rS[i]/length/waveLength;
    rK[i] = amplitudeFactor* rC[i]*(std::sin(arg)     + rB[i] +curv_PID_fac);
    vK[i] = amplitudeFactor*(vC[i]*(std::sin(arg)     + rB[i] +curv_PID_fac)
                            +rC[i]*(std::cos(arg)*darg+ vB[i] +curv_PID_dif));
    assert(not std::isnan(rK[i]));
    assert(not std::isinf(rK[i]));
    assert(not std::isnan(vK[i]));
    assert(not std::isinf(vK[i]));
  }

  //need to come up with something for midline torsion here, use zero values for now
  for(int i=0; i<Nm; ++i)
  {
    //const double s = rS[i];
    //rT[i] = A*(exp(-30*s)-K*exp(-100*(s-0.2)*(s-0.2)));
    rT[i] = 0.0;
    vT[i] = 0.0;
  }

  // solve frenet to compute midline parameters
  //Frenet2D::solve(Nm, rS, rK, vK, rX, rY, vX, vY, norX, norY, vNorX, vNorY);
  Frenet3D::solve(Nm, rS, rK, vK, rT, vT, rX, rY, rZ, vX, vY, vZ, norX, norY, norZ, vNorX, vNorY, vNorZ, binX, binY, binZ, vBinX, vBinY, vBinZ);
}

void StefanFish::save(std::string filename)
{
  auto * const cFish = dynamic_cast<CurvatureDefinedFishData*>( myFish );
  if(cFish == nullptr) { printf("Someone touched my fish\n"); abort(); }

  //assert(std::abs(t-sim_time)<std::numeric_limits<Real>::epsilon());
  std::ofstream savestream;
  savestream.setf(std::ios::scientific);
  savestream.precision(std::numeric_limits<Real>::digits10 + 1);
  savestream.open(filename + ".txt");

  //const double timeshift = myFish->timeshift;
  //const double time0 = myFish->time0;
  //const double l_Tp = myFish->l_Tp;

  savestream<<sim.time<<"\t"<<sim.dt<<std::endl;
  savestream<<position[0]<<"\t"<<position[1]<<"\t"<<position[2]<<std::endl;
  savestream<<quaternion[0]<<"\t"<<quaternion[1]<<"\t"<<quaternion[2]<<"\t"<<quaternion[3]<<std::endl;
  savestream<<transVel[0]<<"\t"<<transVel[1]<<"\t"<<transVel[2]<<std::endl;
  savestream<<angVel[0]<<"\t"<<angVel[1]<<"\t"<<angVel[2]<<std::endl;
  savestream.close();

  cFish->curvatureScheduler.save(filename+"_curv");
  cFish->rlBendingScheduler.save(filename+"_control");
}

void StefanFish::restart(std::string filename)
{
  auto * const cFish = dynamic_cast<CurvatureDefinedFishData*>( myFish );
  if(cFish == nullptr) { printf("Someone touched my fish\n"); abort(); }

  double restarted_time, restarted_dt; // timeshift, time0, l_Tp,
  std::ifstream restartstream;
  restartstream.open(filename+".txt");
  if(!restartstream.good()){
    printf("Could not restart from file\n");
    return;
  }
  restartstream >> restarted_time >> restarted_dt;
  restartstream >> position[0] >> position[1] >> position[2];
  restartstream >> quaternion[0] >> quaternion[1] >> quaternion[2] >> quaternion[3];
  restartstream >> transVel[0] >> transVel[1] >> transVel[2];
  restartstream >> angVel[0] >> angVel[1] >> angVel[2];
  restartstream.close();

  cFish->curvatureScheduler.restart(filename+"_curv");
  cFish->rlBendingScheduler.restart(filename+"_control");

  //myFish->timeshift = timeshift;
  //myFish->time0 = time0;
  //myFish->l_Tp = l_Tp;

  if(!sim.rank)
  {
  std::cout<<"RESTARTED FISH: "<<std::endl;
  std::cout<<"TIME, DT: "<<restarted_time<<" "<<restarted_dt<<std::endl;
  std::cout<<"POS: "<<position[0]<<" "<<position[1]<<" "<<position[2]<<std::endl;
  std::cout<<"ANGLE: "<<quaternion[0]<<" "<<quaternion[1]<<" "<<quaternion[2]<<" "<<quaternion[3]<<std::endl;
  std::cout<<"TVEL: "<<transVel[0]<<" "<<transVel[1]<<" "<<transVel[2]<<std::endl;
  std::cout<<"AVEL: "<<angVel[0]<<" "<<angVel[1]<<" "<<angVel[2]<<std::endl;
  }
}

StefanFish::StefanFish(SimulationData & s, ArgumentParser&p) : Fish(s, p)
{
  const double ampFac = p("-amplitudeFactor").asDouble(1.0);
  myFish = new CurvatureDefinedFishData(length, Tperiod, phaseShift,
    sim.hmin, ampFac);

  std::string heightName = p("-heightProfile").asString("baseline");
  std::string  widthName = p( "-widthProfile").asString("baseline");
  MidlineShapes::computeWidthsHeights(heightName, widthName, length,
    myFish->rS, myFish->height, myFish->width, myFish->Nm, sim.rank);

  origC[0] = position[0];
  origC[1] = position[1];
  origAng = 2 * std::atan2(quaternion[3], quaternion[0]);
  //bool bKillAmplitude = parser("-zeroAmplitude").asInt(0);
  //if(bKillAmplitude) myFish->killAmplitude();

  if(!sim.rank) printf("nMidline=%d, length=%f, Tperiod=%f, phaseShift=%f\n",myFish->Nm, length, Tperiod, phaseShift);
}

//static inline Real sgn(const Real val) { return (0 < val) - (val < 0); }
void StefanFish::create()
{
  auto * const cFish = dynamic_cast<CurvatureDefinedFishData*>( myFish );

  const double DT = sim.dt/Tperiod;
  //const double time = sim.time;
  // Control pos diffs
  const double   xDiff = (position[0] - origC[0])/length;
  const double   yDiff = (position[1] - origC[1])/length;
  const double angDiff =  2 * std::atan2(quaternion[3], quaternion[0])    - origAng;
  const double relU = (transVel[0] + sim.uinf[0]) / length;
  const double relV = (transVel[1] + sim.uinf[1]) / length;
  const double aVelZ = angVel[2], lastAngVel = cFish->lastAvel;
  // compute ang vel at t - 1/2 dt such that we have a better derivative:
  const double aVelMidP = (aVelZ + lastAngVel)*Tperiod/2;
  const double aVelDiff = (aVelZ - lastAngVel)*Tperiod/sim.dt;
  cFish->lastAvel = aVelZ; // store for next time

  // derivatives of following 2 exponential averages:
  const double velDAavg = (angDiff-cFish->avgDangle)/Tperiod + DT * aVelZ;
  const double velDYavg = (  yDiff-cFish->avgDeltaY)/Tperiod + DT * relV;
  const double velAVavg = 10*((aVelMidP-cFish->avgAngVel)/Tperiod +DT*aVelDiff);
  // exponential averages
  cFish->avgDangle = (1.0 -DT) * cFish->avgDangle +    DT * angDiff;
  cFish->avgDeltaY = (1.0 -DT) * cFish->avgDeltaY +    DT *   yDiff;
  // faster average:
  cFish->avgAngVel = (1-10*DT) * cFish->avgAngVel + 10*DT *aVelMidP;
  const double avgDangle = cFish->avgDangle, avgDeltaY = cFish->avgDeltaY;

  // integral (averaged) and proportional absolute DY and their derivative
  const double absPy = std::fabs(yDiff), absIy = std::fabs(avgDeltaY);
  const double velAbsPy =     yDiff>0 ? relV     : -relV;
  const double velAbsIy = avgDeltaY>0 ? velDYavg : -velDYavg;

  if (bCorrectPosition || bCorrectTrajectory)
    assert(origAng<2e-16 && "TODO: rotate pos and vel to fish POV to enable \
                             PID to work even for non-zero angles");

  if (bCorrectPosition && sim.dt>0)
  {
    //If angle is positive: positive curvature only if Dy<0 (must go up)
    //If angle is negative: negative curvature only if Dy>0 (must go down)
    const double IangPdy = (avgDangle *     yDiff < 0)? avgDangle * absPy : 0;
    const double PangIdy = (angDiff   * avgDeltaY < 0)? angDiff   * absIy : 0;
    const double IangIdy = (avgDangle * avgDeltaY < 0)? avgDangle * absIy : 0;

    // derivatives multiplied by 0 when term is inactive later:
    const double velIangPdy = velAbsPy * avgDangle + absPy * velDAavg;
    const double velPangIdy = velAbsIy * angDiff   + absIy * aVelZ;
    const double velIangIdy = velAbsIy * avgDangle + absIy * velDAavg;

    //zero also the derivatives when appropriate
    const double coefIangPdy = avgDangle *     yDiff < 0 ? 1 : 0;
    const double coefPangIdy = angDiff   * avgDeltaY < 0 ? 1 : 0;
    const double coefIangIdy = avgDangle * avgDeltaY < 0 ? 1 : 0;

    const double valIangPdy = coefIangPdy *    IangPdy;
    const double difIangPdy = coefIangPdy * velIangPdy;
    const double valPangIdy = coefPangIdy *    PangIdy;
    const double difPangIdy = coefPangIdy * velPangIdy;
    const double valIangIdy = coefIangIdy *    IangIdy;
    const double difIangIdy = coefIangIdy * velIangIdy;
    const double periodFac = 1.0 - xDiff;
    const double periodVel =     - relU;

    const double totalTerm = valIangPdy + valPangIdy + valIangIdy;
    const double totalDiff = difIangPdy + difPangIdy + difIangIdy;
    cFish->correctTrajectory(totalTerm, totalDiff, sim.time, sim.dt);
    cFish->correctTailPeriod(periodFac, periodVel, sim.time, sim.dt);
  }
  // if absIy<EPS then we have just one fish that the simulation box follows
  // therefore we control the average angle but not the Y disp (which is 0)
  else if (bCorrectTrajectory && sim.dt>0)
  {
    const double avgAngVel = cFish->avgAngVel, absAngVel = std::fabs(avgAngVel);
    const double absAvelDiff = avgAngVel>0? velAVavg : -velAVavg;
    const Real coefInst = angDiff*avgAngVel>0 ? 0.01 : 1, coefAvg = 0.1;
    const Real termInst = angDiff*absAngVel;
    const Real diffInst = angDiff*absAvelDiff + aVelZ*absAngVel;
    const double totalTerm = coefInst*termInst + coefAvg*avgDangle;
    const double totalDiff = coefInst*diffInst + coefAvg*velDAavg;

    cFish->correctTrajectory(totalTerm, totalDiff, sim.time, sim.dt);
  }

  // to debug and check state function, but requires an other obstacle
  //const int indCurrAct = (time + sim.dt)/(Tperiod/2);
  //if(time < indCurrAct*Tperiod/2) state(sim.shapes[0]);

  Fish::create();
}

#if 1
void StefanFish::act(const Real t_rlAction, const std::vector<double>& a) const
{
  auto * const cFish = dynamic_cast<CurvatureDefinedFishData*>( myFish );
  if( cFish == nullptr ) { printf("Someone touched my fish\n"); abort(); }
  cFish->execute(sim.time, t_rlAction, a);
}

double StefanFish::getLearnTPeriod() const
{
  auto * const cFish = dynamic_cast<CurvatureDefinedFishData*>( myFish );
  if( cFish == nullptr ) { printf("Someone touched my fish\n"); abort(); }
  return cFish->periodPIDval;
}

double StefanFish::getPhase(const double t) const
{
  auto * const cFish = dynamic_cast<CurvatureDefinedFishData*>( myFish );
  if( cFish == nullptr ) { printf("Someone touched my fish\n"); abort(); }
  const double T0 = cFish->time0;
  const double Ts = cFish->timeshift;
  const double Tp = cFish->periodPIDval;
  const double arg  = 2*M_PI*((t-T0)/Tp +Ts) + M_PI*phaseShift;
  const double phase = std::fmod(arg, 2*M_PI);
  return (phase<0) ? 2*M_PI + phase : phase;
}

std::vector<double> StefanFish::state() const
{
  auto * const cFish = dynamic_cast<CurvatureDefinedFishData*>( myFish );
  if( cFish == nullptr ) { printf("Someone touched my fish\n"); abort(); }
  std::vector<double> S(10,0);
  S[0] = ( position[0] - origC[0] )/ length;
  S[1] = ( position[1] - origC[1] )/ length;
  S[2] = 2 * std::atan2(quaternion[3], quaternion[0]);
  S[3] = getPhase( sim.time );
  S[4] = transVel[0] * Tperiod / length;
  S[5] = transVel[1] * Tperiod / length;
  S[6] = angVel[2] * Tperiod;
  S[7] = cFish->lastTact;
  S[8] = cFish->lastCurv;
  S[9] = cFish->oldrCurv;

  #ifndef STEFANS_SENSORS_STATE
    return S;
  #else
   S.resize(16); // ONLY 2D shear for consistency with 2D swimmers

    // Get velInfo
    const std::vector<cubism::BlockInfo>& velInfo = sim.vInfo();

    //// Sensor Signal on Front of Fish ////
    ////////////////////////////////////////

    // get front-point
    const std::array<Real,3> pFront = {cFish->sensorLocation[0], cFish->sensorLocation[1], cFish->sensorLocation[2]};

    // first point of the two skins is the same
    // normal should be almost the same: take the mean
    const std::array<Real,3> normalFront = { cFish->sensorNormals[0]
                                             cFish->sensorNormals[1]
                                             cFish->sensorNormals[2]};

    // compute shear stress
    std::array<Real,2> tipShear = getShear( pFront, normalFront, velInfo );

    //// Sensor Signal on Side of Fish ////
    ///////////////////////////////////////

    std::array<Real,2> lowShear, topShear;

    // get index for sensors on the side of head
    int iHeadSide = 0;
    for(int i=0; i<myFish->Nm-1; ++i)
      if( myFish->rS[i] <= 0.04*length && myFish->rS[i+1] > 0.04*length )
        iHeadSide = i;
    assert(iHeadSide>0);

    for(int a = 0; a<2; ++a)
    {
      // get point
      const std::array<Real,3> pSide = {cFish->sensorLocation[(a+1)*3+0], cFish->sensorLocation[(a+1)*3+1], cFish->sensorLocation[(a+1)*3+2]};

      // get normal to surface
      const std::array<Real,3> normSide = {cFish->sensorNormals[(a+1)*3+0], cFish->sensorNormals[(a+1)*3+1], cFish->sensorNormals[(a+1)*3+2]};

      // compute shear stress
      std::array<Real,2> sideShear = getShear( pSide, normSide, velInfo );

      //Michalis: THIS ONLY WORKS FOR A FISH ON A PLANE
      // now figure out how to rotate it along the fish skin for consistency:
      //const Real dX = D->xSurf[iHeadSide+1] - D->xSurf[iHeadSide];
      //const Real dY = D->ySurf[iHeadSide+1] - D->ySurf[iHeadSide];
      const Real dX = cFish->sensorDelta[0];
      const Real dY = cFish->sensorDelta[1];
      const Real proj = dX * normSide[0] - dY * normSide[1];
      const Real tangX = proj>0?  normSide[0] : -normSide[0]; // s.t. tang points from head
      const Real tangY = proj>0? -normSide[1] :  normSide[1]; // to tail, normal outward
      (a==0? topShear[0] : lowShear[0]) = sideShear[0] * normSide[0] + sideShear[1] * normSide[1];
      (a==0? topShear[1] : lowShear[1]) = sideShear[0] * tangX + sideShear[1] * tangY;
    }

    // put non-dimensional results into state into state
    S[10] = tipShear[0] * Tperiod / length;
    S[11] = tipShear[1] * Tperiod / length;
    S[12] = lowShear[0] * Tperiod / length;
    S[13] = lowShear[1] * Tperiod / length;
    S[14] = topShear[0] * Tperiod / length;
    S[15] = topShear[1] * Tperiod / length;
    // printf("shear tip:[%f %f] lower side:[%f %f] upper side:[%f %f]\n", S[10],S[11], S[12],S[13], S[14],S[15]);
    // fflush(0);
    return S;
  #endif
}

/* helpers to compute sensor information */
#ifdef STEFANS_SENSORS_STATE
// function that finds block id of block containing pos (x,y)
ssize_t StefanFish::holdingBlockID(const std::array<Real,3> pos, const std::vector<cubism::BlockInfo>& velInfo) const
{
  for(size_t i=0; i<velInfo.size(); ++i)
  {
    // get gridspacing in block
    const Real h = velInfo[i].h;

    // compute lower left corner of block
    std::array<Real,3> MIN = velInfo[i].pos<Real>(0, 0, 0);
    for(int j=0; j<2; ++j)
      MIN[j] -= 0.5 * h; // pos returns cell centers

    // compute top right corner of block
    std::array<Real,3> MAX = velInfo[i].pos<Real>(FluidBlock::sizeX-1, FluidBlock::sizeY-1, FluidBlock::sizeZ-1);
    for(int j=0; j<2; ++j)
      MAX[j] += 0.5 * h; // pos returns cell centers

    // check whether point is inside block
    if( pos[0] > MIN[0] && pos[1] > MIN[1] && pos[2] > MIN[2] && pos[0] <= MAX[0] && pos[1] <= MAX[1] && pos[2] <= MAX[2] )
    {
      // return blockId holding point
      return i;
    }
  }
  // rank does not contain point
  return -1;
};

// function that gives indice of point in block
std::array<int, 3> StefanFish::safeIdInBlock(const std::array<Real,3> pos, const std::array<Real,3> org, const Real invh ) const
{
  const int indx = (int) std::round((pos[0] - org[0])*invh);
  const int indy = (int) std::round((pos[1] - org[1])*invh);
  const int indz = (int) std::round((pos[2] - org[2])*invh);
  const int ix = std::min( std::max(0, indx), FluidBlock::sizeX-1);
  const int iy = std::min( std::max(0, indy), FluidBlock::sizeY-1);
  const int iz = std::min( std::max(0, indz), FluidBlock::sizeZ-1);
  return std::array<int, 3>{{ix, iy, iz}};
};

// returns shear at given surface location
std::array<Real, 2> StefanFish::getShear(const std::array<Real,3> pSurf, const std::array<Real,3> normSurf, const std::vector<cubism::BlockInfo>& velInfo) const
{
  // Buffer to broadcast x-/y-velocities and gridspacing
  double velocityH[3] = {0.0, 0.0, 0.0};

  // 1. Compute surface velocity on surface
  // get blockId of surface
  const ssize_t blockIdSurf = holdingBlockID(pSurf, velInfo);

  // get surface velocity if block containing point found
  if( blockIdSurf >= 0 ) {
    // check whether obstacle block exists
    if(obstacleBlocks[skinBinfo.blockID] == nullptr ){
      printf("[CUP2D, rank %u] obstacleBlocks[%llu] is a nullptr! obstacleBlocks.size()=%lu", sim.rank, skinBinfo.blockID, obstacleBlocks.size());
      fflush(0);
      abort();
    }

    // get block
    const auto& skinBinfo = velInfo[blockIdSurf];

    // get origin of block
    const std::array<Real,3> oBlockSkin = skinBinfo.pos<Real>(0, 0, 0);

    // get gridspacing on this block
    velocityH[2] = velInfo[blockIdSurf].h;

    // get index of point in block
    const std::array<int,3> iSkin = safeIdInBlock(pSurf, oBlockSkin, 1/velocityH[2]);

    // get deformation velocity
    const Real* const udef = obstacleBlocks[skinBinfo.blockID]->udef[iSkin[2]][iSkin[1]][iSkin[0]];

    // compute velocity of skin point
    velocityH[0]=transVel[0]+( angVel[1]*(pSurf[2]-centerOfMass[2])-angVel[2]*(pSurf[1]-centerOfMass[1]))+udef[0];
    velocityH[1]=transVel[1]+(-angVel[2]*(pSurf[0]-centerOfMass[0])+angVel[0]*(pSurf[2]-centerOfMass[2]))+udef[1];
  }

  // Allreduce to Bcast surface velocity
  MPI_Allreduce(MPI_IN_PLACE, velocityH, 3, MPI_DOUBLE, MPI_SUM, sim.app_comm);

  // Assign skin velocities and grid-spacing
  const double uSkin = velocityH[0];
  const double vSkin = velocityH[1];
  const double h     = velocityH[2];
  const double invh = 1/h;

  // Reset buffer to 0
  velocityH[0] = 0.0; velocityH[1] = 0.0; velocityH[2] = 0.0;

  // 2. Compute flow velocity away from surface
  // compute point on lifted surface
  const std::array<Real,3> pLiftedSurf = { pSurf[0] + h * normSurf[0],
                                           pSurf[1] + h * normSurf[1],
                                           pSurf[2] + h * normSurf[2] };

  // get blockId of lifted surface
  const ssize_t blockIdLifted = holdingBlockID(pLiftedSurf, velInfo);

  // get surface velocity if block containing point found
  if( blockIdLifted >= 0 ) {
    // get block
    const auto& liftedBinfo = velInfo[blockIdLifted];

    // get origin of block
    const std::array<Real,3> oBlockLifted = liftedBinfo.pos<Real>(0, 0, 0);

    // get inverse gridspacing in block
    const Real invhLifted = 1/velInfo[blockIdLifted].h;

    // get index for sensor
    const std::array<int,3> iSens = safeIdInBlock(pLiftedSurf, oBlockLifted, invhLifted);

    // get velocity field at point
    const FluidBlock& b = * (const FluidBlock*) liftedBinfo.ptrBlock;
    velocityH[0] = b(iSens[0], iSens[1], iSens[2]).u;
    velocityH[1] = b(iSens[0], iSens[1], iSens[2]).v;
  }

  // Allreduce to Bcast flow velocity
  MPI_Allreduce(MPI_IN_PLACE, velocityH, 3, MPI_DOUBLE, MPI_SUM, sim.app_comm);

  // Assign lifted skin velocities
  const double uLifted = velocityH[0];
  const double vLifted = velocityH[1];

  // return shear
  return std::array<Real, 2>{{(uLifted - uSkin) * invh,
                              (vLifted - vSkin) * invh }};

};
#endif
#endif

CubismUP_3D_NAMESPACE_END
