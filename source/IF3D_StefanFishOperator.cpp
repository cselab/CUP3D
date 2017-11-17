//
//  CubismUP_3D
//
//  Written by Guido Novati ( novatig@ethz.ch ).
//  This file started as an extension of code written by Wim van Rees
//  Copyright (c) 2017 ETHZ. All rights reserved.
//

#include "IF3D_StefanFishOperator.h"
#include "IF3D_FishLibrary.h"

void IF3D_StefanFishOperator::save(const int step_id, const double t, std::string filename)
{
  //assert(std::abs(t-sim_time)<std::numeric_limits<Real>::epsilon());
  std::ofstream savestream;
  savestream.setf(std::ios::scientific);
  savestream.precision(std::numeric_limits<Real>::digits10 + 1);
  savestream.open(filename + ".txt");

  const double timeshift = myFish->timeshift;
  const double time0 = myFish->time0;
  const double l_Tp = myFish->l_Tp;

  savestream<<t<<"\t"<<sim_dt<<std::endl;
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
  sr.save(step_id, filename);
}

void IF3D_StefanFishOperator::restart(const double t, string filename)
{
    double timeshift, time0, l_Tp;
    std::ifstream restartstream;
    restartstream.open(filename+".txt");
    if(!restartstream.good()){
      printf("Could not restart from file\n");
      return;
    }
    restartstream >> sim_time >> sim_dt;
    assert(std::abs(sim_time-t) < std::numeric_limits<Real>::epsilon());
    restartstream >> position[0] >> position[1] >> position[2];
    restartstream >> quaternion[0] >> quaternion[1] >> quaternion[2] >> quaternion[3];
    restartstream >> transVel[0] >> transVel[1] >> transVel[2];
    restartstream >> angVel[0] >> angVel[1] >> angVel[2];
    restartstream >> theta_internal >> angvel_internal >> adjTh;
    restartstream >> timeshift >> time0 >> l_Tp;// >> new_curv >> old_curv >> new_Tp;
    //restartstream >> _2Dangle >> old_curv >> new_curv;
    restartstream.close();

    sr.restart(filename);
    myFish->curvScheduler.restart(filename+"_curv");
    myFish->baseScheduler.restart(filename+"_base");
    myFish->adjustScheduler.restart(filename+"_adj");
    myFish->timeshift = timeshift;
    myFish->time0 = time0;
    myFish->l_Tp = l_Tp;

    if(!rank)
    {
    std::cout<<"RESTARTED FISH: "<<std::endl;
    std::cout<<"TIME, DT: "<<sim_time<<" "<<sim_dt<<std::endl;
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

IF3D_StefanFishOperator::IF3D_StefanFishOperator(FluidGridMPI*g, ArgumentParser&p, const Real*const u) : IF3D_FishOperator(g, p, u)
{
  _parseArguments(p);
  const int Nextension = NEXTDX*NPPEXT;// up to 3dx on each side (to get proper interpolation up to 2dx)
  const double h = vInfo[0].h_gridpoint;
  const double target_Nm = TGTPPB*length/h;
  const double dx_extension = (1./NEXTDX)*h;
  const int Nm = (Nextension+1)*(int)std::ceil(target_Nm/(Nextension+1)) + 1;

  if(!rank) printf("%d %f %f %f %f\n", Nm, length, Tperiod, phaseShift, dx_extension);

  myFish = new CurvatureDefinedFishData(Nm, length, Tperiod, phaseShift, dx_extension);

  sr.updateInstant(position[0], absPos[0], position[1], absPos[1],
                    _2Dangle, transVel[0], transVel[1], angVel[2]);
}

void IF3D_StefanFishOperator::_parseArguments(ArgumentParser & parser)
{
  IF3D_FishOperator::_parseArguments(parser);
  sr = StateReward(length, Tperiod);
  sr.parseArguments(parser);
}

void IF3D_StefanFishOperator::execute(const int iAgent, const double time, const vector<double> act)
{
  const int nActions = act.size();
  const double eps = std::numeric_limits<Real>::epsilon();
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

  #ifndef __RL_TRAINING
  if(!rank) {
    printf("Next action of agent %d at time %g\n", obstacleID, sr.t_next_comm);
    ofstream filedrag;
    filedrag.open(("orders_"+to_string(obstacleID)+".txt").c_str(), ios::app);
    filedrag<<time<<" "<<act[0];
    if(nActions==2) filedrag<<" "<<act[1];
    filedrag<<endl;
    filedrag.close();
  }
  #endif
}
