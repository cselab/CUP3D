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
  Real * const rK;
  Real * const vK;
  Real * const rC;
  Real * const vC;
  Real * const rB;
  Real * const vB;

 public:

  CurvatureDefinedFishData(double L, double T, double phi, double _h, const double _ampFac)
  : FishMidlineData(L, T, phi, _h, _ampFac),
    rK(_alloc(Nm)),vK(_alloc(Nm)), rC(_alloc(Nm)),vC(_alloc(Nm)),
    rB(_alloc(Nm)),vB(_alloc(Nm)) { }

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
  }

  void computeMidline(const double time, const double dt) override;
};

void CurvatureDefinedFishData::execute(const double time, const double l_tnext,
                                       const std::vector<double>& input)
{
  if (input.size()>1) {
    rlBendingScheduler.Turn(input[0], l_tnext);
    //first, shift time to  previous turn node
    timeshift += (l_tnext-time0)/periodPIDval;
    time0 = l_tnext;
    periodPIDval = Tperiod*(1.+input[1]);
    periodPIDdif = 0;
  } else if (input.size()>0) {
    printf("Turning by %g at time %g with period %g.\n", input[0], time, l_tnext);
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

  //printf("%g %g %g %g\n", rB[12], rB[425], rB[838], rB[1238]);
  #if 0
    { // we dump the profile points
      FILE * f = fopen("stefan.dat","a");
      std::array<Real, 6> curv,base;
      curvScheduler.ParameterScheduler<6>::gimmeValues(time, curv);
      baseScheduler.ParameterScheduler<6>::gimmeValues(time, base);
      fprintf(f,"%9.9e %9.9e %9.9e %9.9e %9.9e %9.9e %9.9e %9.9e %9.9e %9.9e %9.9e %9.9e %9.9e %9.9e\n",
          time,curv[0],curv[1],curv[2],curv[3],curv[4],curv[5],base[0],base[1],base[2],base[3],base[4],base[5]);
      fclose(f);
    }
  #endif
  // solve frenet to compute midline parameters
  Frenet2D::solve(Nm, rS, rK, vK, rX, rY, vX, vY, norX, norY, vNorX, vNorY);
  #if 0
    #warning USED MPI COMM WORLD
    {
      int rank;
      MPI_Comm_rank(MPI_COMM_WORLD,&rank);
      if (rank!=0) return;
      FILE * f = fopen("stefan_profile","w");
      for(int i=0;i<Nm;++i)
        fprintf(f,"%d %g %g %g %g %g %g %g %g %g\n",
          i,rS[i],rX[i],rY[i],vX[i],vY[i],
          vNorX[i],vNorY[i],width[i],height[i]);
      fclose(f);
    }
  #endif
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
  savestream<<theta_internal<<"\t"<<angvel_internal<<std::endl;//<<"\t"<<adjTh
  //savestream<<timeshift<<"\t"<<time0<<"\t"<<l_Tp<<std::endl;
  //<<"\t"<<new_curv<<"\t"<<old_curv<<"\t"<<new_Tp<<std::endl;
  //savestream<<_2Dangle<<"\t"<<old_curv<<"\t"<<new_curv<<std::endl;
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
  restartstream >> theta_internal >> angvel_internal; //  >> adjTh
  //restartstream >> timeshift >> time0 >> l_Tp;// >> new_curv >> old_curv >> new_Tp;
  //restartstream >> _2Dangle >> old_curv >> new_curv;
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
  std::cout<<"INTERN: "<<theta_internal<<" "<<angvel_internal<<std::endl;
  //std::cout<<"TIMESHIFT: "<<timeshift<<" "<<time0<<" "<<l_Tp<<std::endl;
  //std::cout<<"ACTIONS: "<<new_curv<<" "<<old_curv<<" "<<new_Tp<<std::endl;
  std::cout<<"2D angle: "<<_2Dangle<<std::endl;
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
  origAng = _2Dangle;
  //bool bKillAmplitude = parser("-zeroAmplitude").asInt(0);
  //if(bKillAmplitude) myFish->killAmplitude();

  if(!sim.rank) printf("%d %f %f %f\n",myFish->Nm, length, Tperiod, phaseShift);
}

//static inline Real sgn(const Real val) { return (0 < val) - (val < 0); }
void StefanFish::create()
{
  auto * const cFish = dynamic_cast<CurvatureDefinedFishData*>( myFish );
  if(cFish == nullptr) { printf("Someone touched my fish\n"); abort(); }
  const double DT = sim.dt/Tperiod;
  //const double time = sim.time;
  // Control pos diffs
  const double   xDiff = (position[0] - origC[0])/length;
  const double   yDiff = (position[1] - origC[1])/length;
  const double angDiff =  _2Dangle    - origAng;
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

    //if(not sim.muteAll) {
    //  std::ofstream filePID;
    //  std::stringstream ssF;
    //  ssF<<"PID_"<<obstacleID<<".dat";
    //  filePID.open(ssF.str().c_str(), std::ios::app);
    //  filePID<<time<<" "<<valIangPdy<<" "<<difIangPdy
    //               <<" "<<valPangIdy<<" "<<difPangIdy
    //               <<" "<<valIangIdy<<" "<<difIangIdy
    //               <<" "<<periodFac <<" "<<periodVel <<"\n";
    //}
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

    //if(not sim.muteAll) {
    //  std::ofstream filePID;
    //  std::stringstream ssF;
    //  ssF<<"PID_"<<obstacleID<<".dat";
    //  filePID.open(ssF.str().c_str(), std::ios::app);
    //  filePID<<time<<" "<<coefInst*termInst<<" "<<coefInst*diffInst
    //               <<" "<<coefAvg*avgDangle<<" "<<coefAvg*velDAavg<<"\n";
    //}
    cFish->correctTrajectory(totalTerm, totalDiff, sim.time, sim.dt);
  }

  // to debug and check state function, but requires an other obstacle
  //const int indCurrAct = (time + sim.dt)/(Tperiod/2);
  //if(time < indCurrAct*Tperiod/2) state(sim.shapes[0]);

  Fish::create();
}

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
  S[2] = _2Dangle;
  S[3] = getPhase( sim.time );
  S[4] = transVel[0] * Tperiod / length;
  S[5] = transVel[1] * Tperiod / length;
  S[6] = angVel[2] * Tperiod;
  S[7] = cFish->lastTact;
  S[8] = cFish->lastCurv;
  S[9] = cFish->oldrCurv;

  #if 1
    return S;
  #else // TODO
    S.resize(16);

    // Get blocks on this rank
    const std::vector<cubism::BlockInfo>& myInfo = sim.grid->getBlocksInfo();

    b(ix,iy).u[0];
    b(ix,iy).s;

    b(ix,iy,iz).u;
    b(ix,iy,iz).v;
    b(ix,iy,iz).w;
    b(ix,iy,iz).p;

    // Get fish skin
    const auto &DU = myFish->upperSkin, &DL = myFish->lowerSkin;

    //// Sensor Signal on Front of Fish ////
    ////////////////////////////////////////

    std::array<Real,2> tipShear;

    // get front-point
    const std::array<Real,2> pFront = {DU.xSurf[0], DU.ySurf[0]};

    // get blockId for this point
    const int blockIdFront = holdingBlockID(pFront, myInfo);

    if( blockIdFront > 0 )
    {
      // get gridspacing and its inverse in block
      const Real hFront = myInfo[blockIdFront].h_gridpoint, invhFront = 1/hFront;

      // first point of the two skins is the same
      // normal should be almost the same: take the mean
      const Real normXfront = (DU.normXSurf[0] + DL.normXSurf[0]) / 2;
      const Real normYfront = (DU.normYSurf[0] + DL.normYSurf[0]) / 2;

      // compute point for sensor
      const std::array<Real,2> pFrontSens = {pFront[0] + hFront * normXfront,
                                        pFront[1] + hFront * normYfront};\
    }

    // get surface velocity and velocity at sensor
    const std::array<Real,2> vSensFront = sensVel(pFrontSens, myInfo), vSkinFront = skinVel(pFront, myInfo);

    tipShear[0] = (vSensFront[0] - vSkinFront[0]) * invhFront;
    tipShear[1] = (vSensFront[1] - vSkinFront[1]) * invhFront;

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
      // distinguish upper and lower skin
      const auto& D = a==0 ? myFish->upperSkin : myFish->lowerSkin;

      // get point
      const std::array<Real,2> pSide = {D.midX[iHeadSide], D.midY[iHeadSide]};

      // get blockId
      const size_t blockIdSide = holdingBlockID(pSide, velInfo);

      // get gridspacing and its inverse in block
      const Real hSide = velInfo[blockIdSide].h_gridpoint, invhSide = 1/hSide;

      // get normal to surface
      const Real normX = D.normXSurf[iHeadSide], normY = D.normYSurf[iHeadSide];

      // compute point at sensor location
      const std::array<Real,2> pSideSens = {pSide[0] + hSide * normX,
                                        pSide[1] + hSide * normY};

      // get surface velocity and velocity at sensor
      const std::array<Real,2> vSens = sensVel(pSideSens, velInfo), vSkin = skinVel(pSide, velInfo);

      // compute shear
      const Real shearX = (vSens[0] - vSkin[0]) * invhSide;
      const Real shearY = (vSens[1] - vSkin[1]) * invhSide;

      // now figure out how to rotate it along the fish skin for consistency:
      const Real dX = D.xSurf[iHeadSide+1] - D.xSurf[iHeadSide];
      const Real dY = D.ySurf[iHeadSide+1] - D.ySurf[iHeadSide];
      const Real proj = dX * normX - dY * normY;
      const Real tangX = proj>0?  normX : -normX; // s.t. tang points from head
      const Real tangY = proj>0? -normY :  normY; // to tail, normal outward
      (a==0? topShear[0] : lowShear[0]) = shearX * normX + shearY * normY;
      (a==0? topShear[1] : lowShear[1]) = shearX * tangX + shearY * tangY;
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

#if 0 // TODO

/* helpers to compute sensor information */

// function that finds block id of block containing pos (x,y)
int StefanFish::holdingBlockID(const std::array<Real,2> pos, const std::vector<cubism::BlockInfo>& velInfo) const
{
  for(int i=0; i<velInfo.size(); ++i)
  {
    // get gridspacing in block
    const Real h = velInfo[i].h_gridpoint;

    // compute lower left corner of block
    std::array<Real,2> MIN = velInfo[i].pos<Real>(0, 0);
    for(int j=0; j<2; ++j)
      MIN[j] -= 0.5 * h; // pos returns cell centers

    // compute top right corner of block
    std::array<Real,2> MAX = velInfo[i].pos<Real>(VectorBlock::sizeX-1, VectorBlock::sizeY-1);
    for(int j=0; j<2; ++j)
      MAX[j] += 0.5 * h; // pos returns cell centers

    // check whether point is inside block
    if( pos[0] >= MIN[0] && pos[1] >= MIN[1] && pos[0] <= MAX[0] && pos[1] <= MAX[1] )
    {
      // select obstacle block
      if(obstacleBlocks[i] != nullptr ){
        return i;
      }
    }
  }
  // -1 is placeholder for not-finding the block on this rank
  return -1;
};

// function that gives indice of point in block
std::array<int, 2> StefanFish::safeIdInBlock(const std::array<Real,2> pos, const std::array<Real,2> org, const Real invh ) const
{
  const int indx = (int) std::round((pos[0] - org[0])*invh);
  const int indy = (int) std::round((pos[1] - org[1])*invh);
  const int ix = std::min( std::max(0, indx), VectorBlock::sizeX-1);
  const int iy = std::min( std::max(0, indy), VectorBlock::sizeY-1);
  return std::array<int, 2>{{ix, iy}};
};

// return fish velocity at a point on the fish skin:
std::array<Real, 2> StefanFish::skinVel(const std::array<Real,2> pSkin, const std::vector<cubism::BlockInfo>& velInfo) const
{
  // get blockId for pSkin
  const size_t blockId = holdingBlockID(pSkin, velInfo);

  // get block
  const auto& skinBinfo = velInfo[blockId];

  // get origin of block
  const std::array<Real,2> oSkin = skinBinfo.pos<Real>(0, 0);

  // get inverse gridspacing on this block
  const Real invh = 1/velInfo[blockId].h_gridpoint;

  // get index of point in block
  const std::array<int,2> iSkin = safeIdInBlock(pSkin, oSkin, invh);

  //printf("skin pos:[%f %f] -> block org:[%f %f] ind:[%d %d]\n",
  //  pSkin[0], pSkin[1], oSkin[0], oSkin[1], iSkin[0], iSkin[1]);

  // get deformation velocity
  const Real* const udef = obstacleBlocks[skinBinfo.blockID]->udef[iSkin[1]][iSkin[0]];

  // compute velocity of skin point
  const Real uSkin = u - omega * (pSkin[1]-centerOfMass[1]) + udef[0];
  const Real vSkin = v + omega * (pSkin[0]-centerOfMass[0]) + udef[1];

  return std::array<Real, 2>{{uSkin, vSkin}};
};

// return flow velocity at point of flow sensor:
std::array<Real, 2> StefanFish::sensVel(const std::array<Real,2> pSens, const std::vector<cubism::BlockInfo>& velInfo) const
{
  // get blockId
  const size_t blockId = holdingBlockID(pSens, velInfo);

  // get block
  const auto& sensBinfo = velInfo[blockId];

  // get origin of block
  const std::array<Real,2> oSens = sensBinfo.pos<Real>(0, 0);

  // get inverse gridspacing in block
  const Real invh = 1/velInfo[blockId].h_gridpoint;

  // get index for sensor
  const std::array<int,2> iSens = safeIdInBlock(pSens, oSens, invh);

  //printf("sensor pos:[%f %f] -> block org:[%f %f] ind:[%d %d]\n",
  //  pSens[0], pSens[1], oSens[0], oSens[1], iSens[0], iSens[1]);

  // get velocity field at point
  const VectorBlock& b = * (const VectorBlock*) sensBinfo.ptrBlock;

  return std::array<Real, 2>{{b(iSens[0], iSens[1]).u[0],
                              b(iSens[0], iSens[1]).u[1]}};
};

#endif

CubismUP_3D_NAMESPACE_END
