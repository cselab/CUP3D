//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch) and Wim van Rees.
//

#include "obstacles/IF3D_StefanFishOperator.h"

#include "obstacles/extra/IF3D_FishLibrary.h"

#include "Cubism/ArgumentParser.h"

#include <cmath>
#include <array>

class CurvatureDefinedFishData : public FishMidlineData
{
 protected:
  Real * const rK;
  Real * const vK;
  Real * const rC;
  Real * const vC;
  Real * const rB;
  Real * const vB;
  Real * const rA;
  Real * const vA;
  double controlFac = -1, valPID = 0;
  double controlVel = 0, velPID = 0;
 public:

  CurvatureDefinedFishData(double L, double T, double phi, double _h, const double _ampFac)
  : FishMidlineData(L, T, phi, _h, _ampFac),
    rK(_alloc(Nm)),vK(_alloc(Nm)), rC(_alloc(Nm)),vC(_alloc(Nm)),
    rB(_alloc(Nm)),vB(_alloc(Nm)), rA(_alloc(Nm)),vA(_alloc(Nm)) { }

  void _correctTrajectory(const double dtheta, const double vtheta, const double time, double dt) override;

  void _correctAmplitude(double dAmp, double vAmp, const double time, const double dt) override;

  void execute(const double time, const double l_tnext, const std::vector<double>& input) override;

  ~CurvatureDefinedFishData()
  {
    _dealloc(rK);
    _dealloc(vK);
    _dealloc(rC);
    _dealloc(vC);
    _dealloc(rB);
    _dealloc(vB);
    _dealloc(rA);
    _dealloc(vA);
  }

  void computeMidline(const double time, const double dt) override;
};

void CurvatureDefinedFishData::_correctTrajectory(const double dtheta, const double vtheta, const double time, double dt)
{
  velPID = vtheta;
  valPID = dtheta;

  dt = std::max(std::numeric_limits<double>::epsilon(),dt);
  std::array<double, 6> tmp_curv;
  tmp_curv.fill(dtheta);
  //adjustScheduler.transition(time,time,time+2*dt,tmp_curv, true);
  adjustScheduler.transition(time, time-2*dt, time+2*dt, tmp_curv, true);

}

void CurvatureDefinedFishData::_correctAmplitude(double dAmp, double vAmp, const double time, const double dt)
{
  assert(dAmp>0 && dAmp<2); //buhu
  if(dAmp<=0) { dAmp=0; vAmp=0; }
  controlFac = dAmp;
  controlVel = vAmp;

  //TODO actually should be cubic spline!
  //const Real rampUp = time<Tperiod ? time/Tperiod : 1;
  //const Real fac = dAmp*rampUp/length; //curvature is 1/length
  //const std::array<Real ,6> curvature_values = {
  // fac*.82014, fac*1.46515, fac*2.57136, fac*3.75425, fac*5.09147, fac*5.70449
  //};
  //curvScheduler.transition(time,time,time+2*dt, curvature_values, true);
  //curvScheduler.transition(time, time-dt, time+dt, curvature_values);
}

void CurvatureDefinedFishData::execute(const double time, const double l_tnext, const std::vector<double>& input)
{
  if (input.size()>1) {
    baseScheduler.Turn(input[0], l_tnext);
    //first, shift time to  previous turn node
    timeshift += (l_tnext-time0)/l_Tp;
    time0 = l_tnext;
    l_Tp = Tperiod*(1.+input[1]);
  } else if (input.size()>0) {
    printf("Turning by %g at time %g with period %g.\n", input[0], time, l_tnext);
    baseScheduler.Turn(input[0], l_tnext);
  }
}

void CurvatureDefinedFishData::computeMidline(const double time, const double dt)
{
  const double _1oL = 1./length;
  const double _1oT = 1./l_Tp;
  const std::array<double ,6> curvature_points = {
      0, .15*length, .4*length, .65*length, .9*length, length
  };
  const std::array<double,7> baseline_points = {-.5,-.25,0,.25,.5,.75,1.};
  const std::array<double ,6> curvature_values = {
      0.82014*_1oL, 1.46515*_1oL, 2.57136*_1oL,
      3.75425*_1oL, 5.09147*_1oL, 5.70449*_1oL
  };
  //const std::array<double ,6> curvature_values = std::array<double, 6>();
  const std::array<double,6> curvature_zeros = std::array<double, 6>();
  curvScheduler.transition(time,0,Tperiod,curvature_zeros,curvature_values);
  //curvScheduler.transition(time,0,Tperiod,curvature_values,curvature_values);

  // query the schedulers for current values
  curvScheduler.gimmeValues(  time,               curvature_points, Nm, rS, rC, vC);
  baseScheduler.gimmeValues(  time, l_Tp, length, baseline_points,   Nm, rS, rB, vB);
  adjustScheduler.gimmeValues(time,               curvature_points, Nm, rS, rA, vA);
  if(controlFac>0) {
    const double _vA = velPID, _rA = valPID;
    // construct the curvature
    for(int i=0; i<Nm; i++) {
      const double darg = 2.*M_PI* _1oT;
      const double arg  = 2.*M_PI*(_1oT*(time-time0) +timeshift -rS[i]*_1oL/waveLength) + M_PI*phaseShift;
      //amplitudeFactor=0.0 implies towed. A bit scared to directly multiply vK with ampFactor, but ok as long as either 0 or 1
      rK[i] =   amplitudeFactor* rC[i]*(std::sin(arg)     +rB[i]+_rA)*controlFac;
      vK[i] =   amplitudeFactor* (vC[i]*(std::sin(arg)    +rB[i]+_rA)*controlFac
        + rC[i]*(std::cos(arg)*darg+vB[i]+_vA)*controlFac
        + rC[i]*(std::sin(arg)     +rB[i]+_rA)*controlVel);
    }
  } else {
    // construct the curvature
    for(int i=0; i<Nm; i++) {
      const double darg = 2.*M_PI* _1oT;
      const double arg  = 2.*M_PI*(_1oT*(time-time0) +timeshift -rS[i]*_1oL/waveLength) + M_PI*phaseShift;
      //amplitudeFactor=0.0 implies towed. A bit scared to directly multiply vK with ampFactor, but ok as long as either 0 or 1
      rK[i] =   amplitudeFactor* rC[i]*(std::sin(arg)      + rB[i] + rA[i]);
      vK[i] =   amplitudeFactor* (vC[i]*(std::sin(arg)     + rB[i] + rA[i])
        + rC[i]*(std::cos(arg)*darg + vB[i] + vA[i]));
    }
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
  IF2D_Frenet2D::solve(Nm, rS, rK, vK, rX, rY, vX, vY, norX, norY, vNorX, vNorY);
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

void IF3D_StefanFishOperator::save(std::string filename)
{
  //assert(std::abs(t-sim_time)<std::numeric_limits<Real>::epsilon());
  std::ofstream savestream;
  savestream.setf(std::ios::scientific);
  savestream.precision(std::numeric_limits<Real>::digits10 + 1);
  savestream.open(filename + ".txt");

  const double timeshift = myFish->timeshift;
  const double time0 = myFish->time0;
  const double l_Tp = myFish->l_Tp;

  savestream<<sim.time<<"\t"<<sim.dt<<std::endl;
  savestream<<position[0]<<"\t"<<position[1]<<"\t"<<position[2]<<std::endl;
  savestream<<quaternion[0]<<"\t"<<quaternion[1]<<"\t"<<quaternion[2]<<"\t"<<quaternion[3]<<std::endl;
  savestream<<transVel[0]<<"\t"<<transVel[1]<<"\t"<<transVel[2]<<std::endl;
  savestream<<angVel[0]<<"\t"<<angVel[1]<<"\t"<<angVel[2]<<std::endl;
  savestream<<theta_internal<<"\t"<<angvel_internal<<"\t"<<adjTh<<std::endl;
  savestream<<timeshift<<"\t"<<time0<<"\t"<<l_Tp<<std::endl;
  //<<"\t"<<new_curv<<"\t"<<old_curv<<"\t"<<new_Tp<<std::endl;
  //savestream<<_2Dangle<<"\t"<<old_curv<<"\t"<<new_curv<<std::endl;
  savestream.close();

  myFish->curvScheduler.save(filename+"_curv");
  myFish->baseScheduler.save(filename+"_base");
  myFish->adjustScheduler.save(filename+"_adj");
  #ifdef RL_LAYER
    sr.save(step_id, filename);
  #endif
}

void IF3D_StefanFishOperator::restart(std::string filename)
{
  double timeshift, time0, l_Tp, restarted_time, restarted_dt;
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
  restartstream >> theta_internal >> angvel_internal >> adjTh;
  restartstream >> timeshift >> time0 >> l_Tp;// >> new_curv >> old_curv >> new_Tp;
  //restartstream >> _2Dangle >> old_curv >> new_curv;
  restartstream.close();

  myFish->curvScheduler.restart(filename+"_curv");
  myFish->baseScheduler.restart(filename+"_base");
  myFish->adjustScheduler.restart(filename+"_adj");
  myFish->timeshift = timeshift;
  myFish->time0 = time0;
  myFish->l_Tp = l_Tp;

  if(!sim.rank)
  {
  std::cout<<"RESTARTED FISH: "<<std::endl;
  std::cout<<"TIME, DT: "<<restarted_time<<" "<<restarted_dt<<std::endl;
  std::cout<<"POS: "<<position[0]<<" "<<position[1]<<" "<<position[2]<<std::endl;
  std::cout<<"ANGLE: "<<quaternion[0]<<" "<<quaternion[1]<<" "<<quaternion[2]<<" "<<quaternion[3]<<std::endl;
  std::cout<<"TVEL: "<<transVel[0]<<" "<<transVel[1]<<" "<<transVel[2]<<std::endl;
  std::cout<<"AVEL: "<<angVel[0]<<" "<<angVel[1]<<" "<<angVel[2]<<std::endl;
  std::cout<<"INTERN: "<<theta_internal<<" "<<angvel_internal<<std::endl;
  std::cout<<"TIMESHIFT: "<<timeshift<<" "<<time0<<" "<<l_Tp<<std::endl;
  //std::cout<<"ACTIONS: "<<new_curv<<" "<<old_curv<<" "<<new_Tp<<std::endl;
  std::cout<<"2D angle: "<<_2Dangle<<std::endl;
  }
}

IF3D_StefanFishOperator::IF3D_StefanFishOperator(SimulationData & s,
  ArgumentParser&p) : IF3D_FishOperator(s, p)
{
  const double ampFac = p("-amplitudeFactor").asDouble(1.0);
  myFish = new CurvatureDefinedFishData(length, Tperiod, phaseShift,
    vInfo[0].h_gridpoint, ampFac);

  std::string heightName = p("-heightProfile").asString("baseline");
  std::string  widthName = p( "-widthProfile").asString("baseline");
  MidlineShapes::computeWidthsHeights(heightName, widthName, length,
    myFish->rS, myFish->height, myFish->width, myFish->Nm, sim.rank);

  //bool bKillAmplitude = parser("-zeroAmplitude").asInt(0);
  //if(bKillAmplitude) myFish->killAmplitude();

  if(!sim.rank) printf("%d %f %f %f\n",myFish->Nm, length, Tperiod, phaseShift);

  #ifdef RL_LAYER
  sr = StateReward(length, Tperiod);
  sr.parseArguments(p);
  sr.updateInstant(position[0], absPos[0], position[1], absPos[1],
                    _2Dangle, transVel[0], transVel[1], angVel[2]);
  #endif
}


#if 0
void IF3D_StefanFishOperator::execute(const int iAgent, const double time, const std::vector<double> act)
{
  const int nActions = act.size();
  constexpr double eps = std::numeric_limits<Real>::epsilon();
  (void)eps;

  assert(std::fabs(std::cos(0.5*_2Dangle)-quaternion[0]) < 1e-6);
  assert(std::fabs(quaternion[1]) < eps);
  assert(std::fabs(quaternion[2]) < eps);
  assert(std::fabs(std::sin(0.5*_2Dangle)-quaternion[3]) < 1e-6);

  myFish->execute(time, sr.t_next_comm, act);
  sr.old_curv = sr.new_curv;
  sr.new_curv = act[0];
  if(nActions>1) {
      sr.new_Tp = act[1];
      sr.t_next_comm += .5*myFish->l_Tp;
  } else if (nActions==1) {
      sr.t_next_comm += .5*myFish->Tperiod;
  }
  sr.resetAverage();

  /*
  if(!sim.rank) {
    printf("Next action of agent %d at time %g\n", obstacleID, sr.t_next_comm);
    ofstream filedrag;
    filedrag.open(("orders_"+to_string(obstacleID)+".txt").c_str(), ios::app);
    filedrag<<time<<" "<<act[0];
    if(nActions==2) filedrag<<" "<<act[1];
    filedrag<<endl;
    filedrag.close();
  }
  */
}
#endif
