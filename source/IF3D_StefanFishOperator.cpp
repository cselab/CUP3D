//
//  CubismUP_3D
//
//  Written by Guido Novati ( novatig@ethz.ch ).
//  This file started as an extension of code written by Wim van Rees
//  Copyright (c) 2017 ETHZ. All rights reserved.
//

#include "IF3D_StefanFishOperator.h"
#include "IF3D_FishLibrary.h"


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
  void computeWidthsHeights() override;
 public:

  CurvatureDefinedFishData(const int Nm, const double length, const double Tperiod, const double phaseShift, const double dx_ext)
  : FishMidlineData(Nm,length,Tperiod,phaseShift,dx_ext),
    rK(_alloc(Nm)),vK(_alloc(Nm)), rC(_alloc(Nm)),vC(_alloc(Nm)),
    rB(_alloc(Nm)),vB(_alloc(Nm)), rA(_alloc(Nm)),vA(_alloc(Nm)) {
      computeWidthsHeights();
    }

  void _correctTrajectory(const double dtheta, const double vtheta, const double time, double dt) override;

  void _correctAmplitude(double dAmp, double vAmp, const double time, const double dt) override;

  void execute(const double time, const double l_tnext, const vector<double>& input) override;

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

  void computeMidline(const double time) override;
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
  if(dAmp<=0) {
    dAmp=0;
    vAmp=0;
  }
  controlFac = dAmp;
  controlVel = vAmp;
  //const Real rampUp = time<Tperiod ? time/Tperiod : 1; //TODO actually should be cubic spline!
  //const Real fac = dAmp*rampUp/length; //curvature is 1/length
  //const std::array<Real ,6> curvature_values = {
  //  fac*0.82014, fac*1.46515, fac*2.57136, fac*3.75425, fac*5.09147, fac*5.70449
  //};
  //curvScheduler.transition(time,time,time+2*dt, curvature_values, true);
  //curvScheduler.transition(time, time-dt, time+dt, curvature_values);
}

void CurvatureDefinedFishData::execute(const double time, const double l_tnext, const vector<double>& input)
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

void CurvatureDefinedFishData::computeMidline(const double time)
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
  const std::array<double,6> curvature_zeros = std::array<double, 6>();
  curvScheduler.transition(time,0,Tperiod,curvature_zeros,curvature_values);

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
      rK[i] =   rC[i]*(std::sin(arg)     +rB[i]+_rA)*controlFac;
      vK[i] =   vC[i]*(std::sin(arg)     +rB[i]+_rA)*controlFac
        + rC[i]*(std::cos(arg)*darg+vB[i]+_vA)*controlFac
        + rC[i]*(std::sin(arg)     +rB[i]+_rA)*controlVel;
    }
  } else {
    // construct the curvature
    for(int i=0; i<Nm; i++) {
      const double darg = 2.*M_PI* _1oT;
      const double arg  = 2.*M_PI*(_1oT*(time-time0) +timeshift -rS[i]*_1oL/waveLength) + M_PI*phaseShift;
      rK[i] =   rC[i]*(std::sin(arg)      + rB[i] + rA[i]);
      vK[i] =   vC[i]*(std::sin(arg)      + rB[i] + rA[i])
        + rC[i]*(std::cos(arg)*darg + vB[i] + vA[i]);
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

void CurvatureDefinedFishData::computeWidthsHeights()
{
  //Master default is good old Wimmattia shape:
  const int nh = 8;
  const double xh[8] = {0, 0, .2*length, .4*length,
    .6*length, .8*length, length, length};
  // Slim Zebrafish
  const double yh[8] = {0, 5.5e-2*length, 6.8e-2*length, 7.6e-2*length,
    6.4e-2*length, 7.2e-3*length, 1.1e-1*length, 0};

  // Large fin
  //const Real yh[8] = {0, 5.5e-2*length, 1.8e-1*length, 2e-1*length,
  //  6.4e-2*length, 2e-3*length, 3.25e-1*length, 0};

  /*printf("TailFinSize = %f, Wavelength = %f\n", finSize, waveLength);
  fflush(NULL);
  const Real yh[8] = {0, 5.5e-2*length, 1.8e-1*length, 2e-1*length,
    6.4e-2*length, 2e-3*length, finSize*length, 0};*/


  /* // Tuna clone
  const int nh = 9;
  /*const Real xh[9] = {0, 0, 0.2*length, .4*length,
    .6*length, .9*length, .96*length, length, length};
  finSize = 0.23;
  const Real yh[9] = {0, 5e-2*length, 1.4e-1*length, 1.5e-1*length,
    1.1e-1*length, .0*length, 0.2*length, finSize*length, 0};
  printf("WARNING, CHANGED TAIL SHAPE TO INCREASE SURF AREA BY 10%\n");*/
  /*const Real xh[9] = {0, 0, 0.2*length, .4*length,
          .6*length, .9*length, .96*length, length, length};
  finSize = 0.23;
  const Real yh[9] = {0, 5e-2*length, 1.4e-1*length, 1.5e-1*length,
          1.1e-1*length, .0*length, 0.14*length, finSize*length, 0};
  printf("WARNING, CHANGED TAIL FIN HEIGHT 10PERC LARGER\n");
  */

  //What is this?
  //const Real xh[9] = {0, 0, 0.2*length, .4*length, .6*length,
  //        .9*length, .96*length, length, length};
  //finSize = 0.2;
  //const Real yh[9] = {0, 5e-2*length, 1.4e-1*length, 1.5e-1*length, //1.1e-1*length, .0*length, 0.1*length, finSize*length, 0};
  //printf("WARNING, CHANGED TAIL FIN HEIGHT EQUAL\n");
  //printf("TailFinSize = %f, Wavelength = %f\n", finSize, waveLength);
  //fflush(NULL);

  //Master default is good old Wimmattia shape:
  const int nw = 6;
  const Real xw[6] = {0, 0, length/3., 2*length/3., length, length};
  //const Real yw[6] = {0, 8.9e-2*length, 7.0e-2*length,
  //  3.0e-2*length, 2.0e-2*length, 0};
  const Real yw[6] = {0, 8.9e-2*length, 1.7e-2*length,
    1.6e-2*length, 1.3e-2*length, 0};

  //const double tNACA = 0.15;
  //printf("NACA profile thickness = %f\n", tNACA);
  //for(int i=0;i<Nm;++i) width[i]  = _naca_width(rS[i],length,tNACA);

  //Actually do the b spline integration:
  integrateBSpline(height, xh, yh, nh);
  integrateBSpline(width,  xw, yw, nw);

  //for(int i=0;i<Nm;++i) {
  //  width[i]  = _width(rS[i],length);
  //  height[i] = _height(rS[i],length);
  //}

  // output these suckers
  FILE * heightWidth;
  heightWidth = fopen("widthHeight.txt","w");
  {
    for(int i=0;i<Nm;++i) {
      fprintf(heightWidth, "%f \t %f \t %f \n", rS[i], width[i], height[i]);
    }
  }
  fclose(heightWidth);
}

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
