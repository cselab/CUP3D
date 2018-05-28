//
//  CubismUP_3D
//
//  Written by Guido Novati ( novatig@ethz.ch ).
//  This file started as an extension of code written by Wim van Rees
//  Copyright (c) 2017 ETHZ. All rights reserved.
//

#include "IF3D_DeadFishOperator.h"
#include "IF3D_FishLibrary.h"

IF3D_DeadFishOperator::IF3D_DeadFishOperator(FluidGridMPI*g, ArgumentParser&p, const Real*const u) : IF3D_FishOperator(g, p, u), ext_pos{0,0,0}
{
  _parseArguments(p);
  const int Nextension = NEXTDX*NPPEXT;// up to 3dx on each side (to get proper interpolation up to 2dx)
  const double target_Nm = TGTPPB*length/vInfo[0].h_gridpoint;
  const double dx_extension = (1./NEXTDX)*vInfo[0].h_gridpoint;
  const int Nm = (Nextension+1)*(int)std::ceil(target_Nm/(Nextension+1)) + 1;

  printf("%d %f %f %f %f\n",Nm,length,Tperiod,phaseShift,dx_extension);
  myFish = new CarlingFishMidlineData(Nm, length, Tperiod, phaseShift, dx_extension, 0.);
}

void IF3D_DeadFishOperator::_parseArguments(ArgumentParser & parser)
{
  IF3D_FishOperator::_parseArguments(parser);
  parser.set_strict_mode();
  ID = parser("-DID").asInt();
  
  bForcedInSimFrame[0] = true;
  bForcedInSimFrame[1] = true;
  bForcedInSimFrame[2] = true;

  parser.unset_strict_mode();
  P0 = parser("-xpos").asDouble(0.0);
  Y0 = parser("-ypos").asDouble(0.0);
  Ltow = parser("-Ltow").asDouble(2.5);
  Ttow = parser("-Ttow").asDouble(1.0);
  Atow = parser("-Atow").asDouble(-1.);

  if (Atow>0) {
    transVel[0] = 2*Ltow*length/(Ttow*Tperiod);
    if (ID==0) position[1] -= Atow*length;
    if (ID==1) position[1] += Atow*length;
    if (ID==1) transVel[0] -= 4*Ltow*length/(Ttow*Tperiod);
  }

  ext_pos[0] = position[0];
  ext_pos[1] = position[1];
  ext_pos[2] = position[2];
  printf("created IF2D_DeadFish: xpos=%f ypos=%f L=%f T=%f\n",
      position[0],position[1],length,Tperiod);
  printf("P0=%3.3f Ltow=%3.3f Ttow=%3.3f Atow=%3.3f ID=%d \n",
      P0,Ltow,Ttow,Atow,ID);
}

void IF3D_DeadFishOperator::update(const int step_id, const double t, const double dt, const Real *Uinf)
{
  if (Atow>0) {
    //constant acceleration
    double accel = 4*Ltow*length/(Ttow*Tperiod)/(Ttow*Tperiod);
    //alternates between negative and positive
    if (fmod(t/Tperiod,2.*Ttow)>Ttow) accel=-accel;
    double s_c = 1.0;
    //2 obstacles in antiphase
    if (ID==0) { accel=-accel; s_c=-s_c; }
    //integrate in time constant accel:
    position[0] += dt*transVel[0] + 0.5*accel*dt*dt;
    transVel[0] += accel*dt;
    //why here?
    ext_pos[0] += dt*(transVel[0]-2*Ltow*length/(Ttow*Tperiod)) + 0.5*accel*dt*dt;
    const double arg = .5*M_PI*(P0-ext_pos[0])/Ltow/length;
    position[1] = Y0 + s_c*Atow*length*std::cos(arg);
    const double fac1 = .5*s_c*Atow*M_PI/Ltow*sin(arg);
    transVel[1] = fac1*transVel[0];
    ext_pos[1] += dt*(transVel[1]-Uinf[1]);
    const double angle = atan(fac1);
    const double fac2  = -.25*Atow*s_c*M_PI*M_PI/Ltow*cos(arg)*transVel[0]/Ltow/length;
    angVel[2] = fac2/(1+fac1*fac1);
    /*
    const Real dqdt[4] = {
        0.5*(-angVel[2]*quaternion[3]),
        0.5*(-angVel[2]*quaternion[2]),
        0.5*(+angVel[2]*quaternion[1]),
        0.5*(+angVel[2]*quaternion[0])
    };
    const Real deltaq[4] = {dqdt[0]*dt, dqdt[1]*dt, dqdt[2]*dt, dqdt[3]*dt};
    const Real deltaq_length = std::sqrt(deltaq[0]*deltaq[0]+deltaq[1]*deltaq[1]+deltaq[2]*deltaq[2]+deltaq[3]*deltaq[3]);
    if(deltaq_length>std::numeric_limits<Real>::epsilon()) {
      const Real tanfac = std::tan(deltaq_length)/deltaq_length;
      const Real num[4] = {
          quaternion[0]+tanfac*deltaq[0],
          quaternion[1]+tanfac*deltaq[1],
          quaternion[2]+tanfac*deltaq[2],
          quaternion[3]+tanfac*deltaq[3],
      };

      const Real invDenum = 1./(std::sqrt(num[0]*num[0]+num[1]*num[1]+num[2]*num[2]+num[3]*num[3]));
      quaternion[0] = num[0]*invDenum;
      quaternion[1] = num[1]*invDenum;
      quaternion[2] = num[2]*invDenum;
      quaternion[3] = num[3]*invDenum;
    }
     */
    quaternion[0] = std::cos(0.5*angle);
    quaternion[1] = 0;
    quaternion[2] = 0;
    quaternion[3] = std::sin(0.5*angle);
    printf("accel %f, posx %f, posy %f. velx %f, vely %f, angvel %f\n",
        accel,position[0],position[1],transVel[0],transVel[1],angVel[2]);
    //position[1] = 0.5 + Yamplit*sin(2*M_PI*(position[0]-P0)/Yperiod + M_PI*Yphase);
    //angle = (bTilt) ? angle + atan(2*M_PI*Yamplit*cos(2*M_PI*(position[0]-P0)/Yperiod + M_PI*Yphase)/Yperiod) : angle;
    //cout << Atow << length << P0-position[0]<<endl;

    _writeComputedVelToFile(step_id, t, Uinf);
    if (position[0]<0.1) abort();
  }
}

void IF3D_DeadFishOperator::save(const int stepID, const double t, string filename)
{
    //assert(std::abs(t-sim_time)<std::numeric_limits<Real>::epsilon());
    std::ofstream savestream;
    savestream.setf(std::ios::scientific);
    savestream.precision(std::numeric_limits<Real>::digits10 + 1);
    savestream.open(filename + ".txt");

    savestream<<t<<"\t"<<sim_dt<<std::endl;
    savestream<<position[0]<<"\t"<<position[1]<<"\t"<<position[2]<<std::endl;
    savestream<<ext_pos[0]<<"\t"<<ext_pos[1]<<"\t"<<ext_pos[2]<<std::endl;
    savestream<<quaternion[0]<<"\t"<<quaternion[1]<<"\t"<<quaternion[2]<<"\t"<<quaternion[3]<<std::endl;
    savestream<<transVel[0]<<"\t"<<transVel[1]<<"\t"<<transVel[2]<<std::endl;
    savestream<<angVel[0]<<"\t"<<angVel[1]<<"\t"<<angVel[2]<<std::endl;
    savestream<<theta_internal<<"\t"<<angvel_internal<<"\t"<<adjTh<<std::endl;
    savestream.close();

}

void IF3D_DeadFishOperator::restart(const double t, string filename)
{
    std::ifstream restartstream;
    restartstream.open(filename+".txt");
    restartstream >> sim_time >> sim_dt;
    assert(std::abs(sim_time-t) < std::numeric_limits<Real>::epsilon());
    restartstream >> position[0] >> position[1] >> position[2];
    restartstream >> ext_pos[0] >> ext_pos[1] >> ext_pos[2];
    restartstream >> quaternion[0] >> quaternion[1] >> quaternion[2] >> quaternion[3];
    restartstream >> transVel[0] >> transVel[1] >> transVel[2];
    restartstream >> angVel[0] >> angVel[1] >> angVel[2];
    restartstream >> theta_internal >> angvel_internal >> adjTh;
    restartstream.close();

  std::cout<<"RESTARTED FISH: "<<std::endl;
  std::cout<<"TIME, DT: "<<sim_time<<" "<<sim_dt<<std::endl;
  std::cout<<"POS: "<<position[0]<<" "<<position[1]<<" "<<position[2]<<std::endl;
  std::cout<<"ANGLE: "<<quaternion[0]<<" "<<quaternion[1]<<" "<<quaternion[2]<<" "<<quaternion[3]<<std::endl;
  std::cout<<"TVEL: "<<transVel[0]<<" "<<transVel[1]<<" "<<transVel[2]<<std::endl;
  std::cout<<"AVEL: "<<angVel[0]<<" "<<angVel[1]<<" "<<angVel[2]<<std::endl;
  std::cout<<"INTERN: "<<theta_internal<<" "<<angvel_internal<<std::endl;
}
