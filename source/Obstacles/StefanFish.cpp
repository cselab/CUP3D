//
//  Cubism3D
//  Copyright (c) 2022 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//

#include "StefanFish.h"
#include "FishLibrary.h"
#include "FishShapes.h"

#include <Cubism/ArgumentParser.h>

#include <array>
#include <cmath>
#include <sstream>
#include <iomanip>

CubismUP_3D_NAMESPACE_BEGIN
using namespace cubism;

void CurvatureDefinedFishData::execute(const Real time, const Real l_tnext, const std::vector<Real>& input)
{
  if (input.size() == 1)
  {
    //1.midline curvature
    rlBendingScheduler.Turn(input[0], l_tnext);
  }
  else if (input.size() == 3)
  {
    assert (control_torsion == false);

    //1.midline curvature
    rlBendingScheduler.Turn(input[0], l_tnext);

    //2.swimming period
    if (TperiodPID) std::cout << "Warning: PID controller should not be used with RL." << std::endl;
    current_period = periodPIDval;
    next_period = Tperiod * (1 + input[1]);
    transition_start = l_tnext;

    //3.pitching motion
    if (time > Tman_finish)
    {
      Tman_start = time;
      Tman_finish = time + 0.25*Tperiod;
      Lman = 0;
      if (std::fabs(input[2]) > 0.01) Lman = 1.0/input[2];
    }
  }
  else if (input.size() == 5)
  {
    assert (control_torsion == true);

    //1.midline curvature
    rlBendingScheduler.Turn(input[0], l_tnext);

    //2.swimming period
    if (TperiodPID) std::cout << "Warning: PID controller should not be used with RL." << std::endl;
    current_period = periodPIDval;
    next_period = Tperiod * (1 + input[1]);
    transition_start = l_tnext;

    //3.midline torsion
    for (int i = 0 ; i < 3 ; i++)
    {
	    torsionValues_previous [i] = torsionValues[i];
	    torsionValues[i] = input[i+2];
    }
    Ttorsion_start = time;
  }
}

void CurvatureDefinedFishData::computeMidline(const Real t, const Real dt)
{
  periodScheduler.transition(t,transition_start,transition_start+transition_duration,current_period,next_period);
  periodScheduler.gimmeValues(t,periodPIDval,periodPIDdif);
  if (transition_start < t && t < transition_start+transition_duration)//timeshift also rampedup
  {
    timeshift = (t - time0)/periodPIDval + timeshift;
    time0 = t;
  }

  const std::array<Real,6> curvaturePoints = { (Real)0         , (Real).15*length, (Real).4*length, 
                                               (Real).65*length, (Real).9 *length, length};
  const std::array<Real,7>      bendPoints = { (Real)-.5, (Real)-.25, (Real)0, (Real).25, 
                                               (Real).5 , (Real) .75, (Real)1};
  const std::array<Real,6> curvatureValues = {
      (Real)0.82014/length, (Real)1.46515/length, (Real)2.57136/length,
      (Real)3.75425/length, (Real)5.09147/length, (Real)5.70449/length};

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

  #pragma omp parallel for
  for(int i=0; i<Nm; ++i)
  {
    const Real arg = arg0 - 2*M_PI*rS[i]/length/waveLength;
    rK[i] = amplitudeFactor* rC[i]*(std::sin(arg)     + rB[i] +curv_PID_fac);
    vK[i] = amplitudeFactor*(vC[i]*(std::sin(arg)     + rB[i] +curv_PID_fac)
                            +rC[i]*(std::cos(arg)*darg+ vB[i] +curv_PID_dif));
    assert(!std::isnan(rK[i]));
    assert(!std::isinf(rK[i]));
    assert(!std::isnan(vK[i]));
    assert(!std::isinf(vK[i]));
  }

  //could come up with something for midline torsion here, use zero values for now
  for(int i=0; i<Nm; ++i)
  {
    rT[i] = 0;
    vT[i] = 0;
  }
  if (control_torsion)
  {
    const std::array<Real ,3> torsionPoints = { 0.0, 0.5*length, length };
    torsionScheduler.transition (t,Ttorsion_start,Ttorsion_start+0.10*Tperiod,torsionValues_previous,torsionValues);
    torsionScheduler.gimmeValues(t,torsionPoints,Nm,rS,rT,vT);
  }

  // solve frenet to compute midline parameters
  //Frenet2D::solve(Nm, rS, rK, vK, rX, rY, vX, vY, norX, norY, vNorX, vNorY);
  Frenet3D::solve(Nm, rS, rK, vK, rT, vT, rX, rY, rZ, vX, vY, vZ, norX, norY, norZ, vNorX, vNorY, vNorZ, binX, binY, binZ, vBinX, vBinY, vBinZ);

  if (std::fabs(Lman) > 1e-6*length && control_torsion == false)
    performPitchingMotion(t);
}

void CurvatureDefinedFishData::performPitchingMotion(const Real t)
{
  const int sgn = Lman > 0 ? +1:-1;
  const Real Rbig = sgn * 1000.0 * length; //"infinite" radius for the cylinder 
  const Real p    = 0.2;//For the first 20% of the total pitching motion time, the cylinder radius 
                        //goes from Rbig to Lman*length. For the next 60% it stays constant and for
                        //the last 20% it transitions back to Rbig. 
  if (t > Tman_finish || t < Tman_start) {Lman = 0; return;}

  //The Scheduler class only works with std::array, so we had to do the following:
  std::array<Real,1> R;
  std::array<Real,1> Rdot;
  const Real dT = Tman_finish - Tman_start;

  if ( t <= Tman_start + p * dT ) 
  {
     const std::array<Real,1> RZero {Rbig};
     const std::array<Real,1> Rvalue{Lman*length};
     turnZScheduler.transition(t, Tman_start, Tman_start + p*dT, RZero, Rvalue);
  }
  else if ( t >= Tman_finish - p*dT)
  {
     const std::array<Real,1> Rvalue{Rbig};
     const std::array<Real,1> RZero{Lman*length};
     turnZScheduler.transition(t,Tman_finish - p*dT,Tman_finish, RZero, Rvalue);
  }

  turnZScheduler.gimmeValues(t,R,Rdot); 

  //for a given cylinder with radius R and Rdot = dR/dt, compute the new midline
  //we do this by wrapping it around the cylinder. 
  //All y-related quantities do not change. The last point (tail of fish, s=Nm-1) does not move.
  for(int i=0; i<Nm; i++)
  {
    const Real d     = length - rS[i];
    const Real theta = d/R[0];
    const Real cos_theta = cos(theta);
    const Real sin_theta = sin(theta);

    rX[i] = rX[i] + d-R[0]*sin_theta;
    rZ[i] = rZ[i] + R[0]*(1.0-cos_theta);

    //This is how a vector (vX,0) gets rotated and wrapped around the cylinder
    const Real vxF =  vX[i]*cos_theta;
    const Real vzF = -vX[i]*sin_theta;

    //The final velocities are the rotated vector plus the d/dt contribution of R(t)
    vX[i] = vxF + Rdot[0]*(theta*cos_theta - sin_theta);
    vZ[i] = vzF + Rdot[0]*(1.0 - cos_theta - sin_theta * theta);
  }

  recomputeNormalVectors();
}

void CurvatureDefinedFishData::recomputeNormalVectors()
{
  //compute normal and binormal vectors for a given midline
  for(int i=1; i<Nm-1; i++)
  {
    //2nd order finite difference for non-uniform grid
    const Real hp = rS[i+1]-rS[i];
    const Real hm = rS[i]-rS[i-1];
    const Real frac = hp/hm;
    const Real am = -frac*frac;
    const Real a  =  frac*frac -1.0;
    const Real ap = 1.0;
    const Real denom = 1.0 / ( hp * (1.0 + frac) );
    const Real tX  = (am * rX[i-1] + a * rX[i] + ap * rX[i+1])*denom;
    const Real tY  = (am * rY[i-1] + a * rY[i] + ap * rY[i+1])*denom;
    const Real tZ  = (am * rZ[i-1] + a * rZ[i] + ap * rZ[i+1])*denom;
    const Real dtX = (am * vX[i-1] + a * vX[i] + ap * vX[i+1])*denom;
    const Real dtY = (am * vY[i-1] + a * vY[i] + ap * vY[i+1])*denom;
    const Real dtZ = (am * vZ[i-1] + a * vZ[i] + ap * vZ[i+1])*denom;
    const Real BDx = norX[i];
    const Real BDy = norY[i];
    const Real BDz = norZ[i];
    const Real dBDx = vNorX[i];
    const Real dBDy = vNorY[i];
    const Real dBDz = vNorZ[i];
    const Real dot = BDx*tX+BDy*tY+BDz*tZ;
    const Real ddot = dBDx*tX+dBDy*tY+dBDz*tZ + BDx*dtX+BDy*dtY+BDz*dtZ;

    //Project the normal vector computed by the Frenet equations onto (-ty,tx,tz) which is
    //a vector on the plane that is perpendicular to the tangent vector t.
    //This projection defines the new normal vector.
    norX[i] = BDx - dot*tX;
    norY[i] = BDy - dot*tY;
    norZ[i] = BDz - dot*tZ;
    const Real inormn = 1.0/sqrt(norX[i]*norX[i]+norY[i]*norY[i]+norZ[i]*norZ[i]);
    norX[i] *= inormn;
    norY[i] *= inormn;
    norZ[i] *= inormn;
    vNorX[i] = dBDx - ddot*tX - dot*dtX;
    vNorY[i] = dBDy - ddot*tY - dot*dtY;
    vNorZ[i] = dBDz - ddot*tZ - dot*dtZ;
    //Compute the bi-normal vector as t x n
    binX[i] =  tY*norZ[i]-tZ*norY[i];
    binY[i] =  tZ*norX[i]-tX*norZ[i];
    binZ[i] =  tX*norY[i]-tY*norX[i];
    const Real inormb = 1.0/sqrt(binX[i]*binX[i]+binY[i]*binY[i]+binZ[i]*binZ[i]);
    binX[i] *= inormb;
    binY[i] *= inormb;
    binZ[i] *= inormb;
    vBinX[i] =  (dtY*norZ[i]+tY*vNorZ[i])-(dtZ*norY[i]+tZ*vNorY[i]);
    vBinY[i] =  (dtZ*norX[i]+tZ*vNorX[i])-(dtX*norZ[i]+tX*vNorZ[i]);
    vBinZ[i] =  (dtX*norY[i]+tX*vNorY[i])-(dtY*norX[i]+tY*vNorX[i]);
  }

  //take care of first and last point
  for (int i = 0 ; i <= Nm-1; i+= Nm-1)
  {
      const int ipm = (i == Nm - 1) ? i-1:i+1;
      const Real ids  = 1.0/(rS[ipm]-rS[i]);
      const Real tX  = (rX[ipm]-rX[i])*ids;
      const Real tY  = (rY[ipm]-rY[i])*ids;
      const Real tZ  = (rZ[ipm]-rZ[i])*ids;
      const Real dtX = (vX[ipm]-vX[i])*ids;
      const Real dtY = (vY[ipm]-vY[i])*ids;
      const Real dtZ = (vZ[ipm]-vZ[i])*ids;

      const Real BDx = norX[i];
      const Real BDy = norY[i];
      const Real BDz = norZ[i];
      const Real dBDx = vNorX[i];
      const Real dBDy = vNorY[i];
      const Real dBDz = vNorZ[i];
      const Real dot = BDx*tX+BDy*tY+BDz*tZ;
      const Real ddot = dBDx*tX+dBDy*tY+dBDz*tZ + BDx*dtX+BDy*dtY+BDz*dtZ;

      norX[i] = BDx - dot*tX;
      norY[i] = BDy - dot*tY;
      norZ[i] = BDz - dot*tZ;
      const Real inormn = 1.0/sqrt(norX[i]*norX[i]+norY[i]*norY[i]+norZ[i]*norZ[i]);
      norX[i] *= inormn;
      norY[i] *= inormn;
      norZ[i] *= inormn;
      vNorX[i] = dBDx - ddot*tX - dot*dtX;
      vNorY[i] = dBDy - ddot*tY - dot*dtY;
      vNorZ[i] = dBDz - ddot*tZ - dot*dtZ;

      binX[i] =  tY*norZ[i]-tZ*norY[i];
      binY[i] =  tZ*norX[i]-tX*norZ[i];
      binZ[i] =  tX*norY[i]-tY*norX[i];
      const Real inormb = 1.0/sqrt(binX[i]*binX[i]+binY[i]*binY[i]+binZ[i]*binZ[i]);
      binX[i] *= inormb;
      binY[i] *= inormb;
      binZ[i] *= inormb;
      vBinX[i] =  (dtY*norZ[i]+tY*vNorZ[i])-(dtZ*norY[i]+tZ*vNorY[i]);
      vBinY[i] =  (dtZ*norX[i]+tZ*vNorX[i])-(dtX*norZ[i]+tX*vNorZ[i]);
      vBinZ[i] =  (dtX*norY[i]+tX*vNorY[i])-(dtY*norX[i]+tY*vNorX[i]);
  }
}

void StefanFish::saveRestart( FILE * f )
{
  assert(f != NULL);
  Fish::saveRestart(f);
  CurvatureDefinedFishData * const cFish = dynamic_cast<CurvatureDefinedFishData*>( myFish );
  std::stringstream ss;
  ss<<std::setfill('0')<<std::setw(7)<<"_"<<obstacleID<<"_";
  std::string filename = "Schedulers"+ ss.str() + ".restart";
  {
     std::ofstream savestream;
     savestream.setf(std::ios::scientific);
     savestream.precision(std::numeric_limits<Real>::digits10 + 1);
     savestream.open(filename);
     {
       const auto & c = cFish->curvatureScheduler;
       savestream << c.t0 << "\t" << c.t1 << std::endl;
       for(int i=0;i<c.npoints;++i)
         savestream << c.parameters_t0[i]  << "\t"
                    << c.parameters_t1[i]  << "\t"
                    << c.dparameters_t0[i] << std::endl;
     }
     {
       const auto & c = cFish->periodScheduler;
       savestream << c.t0 << "\t" << c.t1 << std::endl;
       for(int i=0;i<c.npoints;++i)
         savestream << c.parameters_t0[i]  << "\t"
                    << c.parameters_t1[i]  << "\t"
                    << c.dparameters_t0[i] << std::endl;
     }
     {
      const auto & c = cFish->rlBendingScheduler;
       savestream << c.t0 << "\t" << c.t1 << std::endl;
       for(int i=0;i<c.npoints;++i)
         savestream << c.parameters_t0[i]  << "\t"
                    << c.parameters_t1[i]  << "\t"
                    << c.dparameters_t0[i] << std::endl;
     }
     {
       const auto & c = cFish->torsionScheduler;
       savestream << c.t0 << "\t" << c.t1 << std::endl;
       for(int i=0;i<c.npoints;++i)
         savestream << c.parameters_t0[i]  << "\t"
                    << c.parameters_t1[i]  << "\t"
                    << c.dparameters_t0[i] << std::endl;
     }
     {
       const auto & c = cFish->turnZScheduler;
       savestream << c.t0 << "\t" << c.t1 << std::endl;
       for(int i=0;i<c.npoints;++i)
         savestream << c.parameters_t0[i]  << "\t"
                    << c.parameters_t1[i]  << "\t"
                    << c.dparameters_t0[i] << std::endl;
     }
     savestream.close();
  }

  //Save these numbers for PID controller and other stuff. Maybe not all of them are needed
  //but we don't care, it's only a few numbers.
  fprintf(f, "origC_x: %20.20e\n", (double)origC[0]);
  fprintf(f, "origC_y: %20.20e\n", (double)origC[1]);
  fprintf(f, "origC_z: %20.20e\n", (double)origC[2]);
  fprintf(f, "errP: %20.20e\n", (double)cFish->errP);
  fprintf(f, "errD: %20.20e\n", (double)cFish->errD);
  fprintf(f,"curv_PID_fac             : %20.20e\n",(double)cFish->curv_PID_fac             );
  fprintf(f,"curv_PID_dif             : %20.20e\n",(double)cFish->curv_PID_dif             );
  fprintf(f,"avgDeltaY                : %20.20e\n",(double)cFish->avgDeltaY                );
  fprintf(f,"avgDangle                : %20.20e\n",(double)cFish->avgDangle                );
  fprintf(f,"avgAngVel                : %20.20e\n",(double)cFish->avgAngVel                );
  fprintf(f,"lastTact                 : %20.20e\n",(double)cFish->lastTact                 );
  fprintf(f,"lastCurv                 : %20.20e\n",(double)cFish->lastCurv                 );
  fprintf(f,"oldrCurv                 : %20.20e\n",(double)cFish->oldrCurv                 );
  fprintf(f,"periodPIDval             : %20.20e\n",(double)cFish->periodPIDval             );
  fprintf(f,"periodPIDdif             : %20.20e\n",(double)cFish->periodPIDdif             );
  fprintf(f,"time0                    : %20.20e\n",(double)cFish->time0                    );
  fprintf(f,"timeshift                : %20.20e\n",(double)cFish->timeshift                );
  fprintf(f,"lastTime                 : %20.20e\n",(double)cFish->lastTime                 );
  fprintf(f,"lastAvel                 : %20.20e\n",(double)cFish->lastAvel                 );
  fprintf(f,"Ttorsion_start           : %20.20e\n",(double)cFish->Ttorsion_start           );
  fprintf(f,"Tman_start               : %20.20e\n",(double)cFish->Tman_start               );
  fprintf(f,"Tman_finish              : %20.20e\n",(double)cFish->Tman_finish              );
  fprintf(f,"Lman;                    : %20.20e\n",(double)cFish->Lman                     );
  fprintf(f,"current_period           : %20.20e\n",(double)cFish->current_period           );
  fprintf(f,"next_period              : %20.20e\n",(double)cFish->next_period              );
  fprintf(f,"transition_start         : %20.20e\n",(double)cFish->transition_start         );
  fprintf(f,"transition_duration      : %20.20e\n",(double)cFish->transition_duration      );
  fprintf(f,"torsionValues[0]         : %20.20e\n",(double)cFish->torsionValues[0]         );
  fprintf(f,"torsionValues[1]         : %20.20e\n",(double)cFish->torsionValues[1]         );
  fprintf(f,"torsionValues[2]         : %20.20e\n",(double)cFish->torsionValues[2]         );
  fprintf(f,"torsionValues_previous[0]: %20.20e\n",(double)cFish->torsionValues_previous[0]);
  fprintf(f,"torsionValues_previous[1]: %20.20e\n",(double)cFish->torsionValues_previous[1]);
  fprintf(f,"torsionValues_previous[2]: %20.20e\n",(double)cFish->torsionValues_previous[2]);
  fprintf(f,"TperiodPID               : %d\n"     ,(int)   cFish->TperiodPID               );
  fprintf(f,"control_torsion          : %d\n"     ,(int)   cFish->control_torsion          );
}

void StefanFish::loadRestart( FILE * f )
{
  assert(f != NULL);
  Fish::loadRestart(f);
  CurvatureDefinedFishData* const cFish = dynamic_cast<CurvatureDefinedFishData*>( myFish );
  std::stringstream ss;
  ss<<std::setfill('0')<<std::setw(7)<<"_"<<obstacleID<<"_";
  std::ifstream restartstream;
  std::string filename = "Schedulers"+ ss.str() + ".restart";
  restartstream.open(filename);
  {
     auto & c = cFish->curvatureScheduler;
     restartstream >> c.t0 >> c.t1;
     for(int i=0;i<c.npoints;++i)
       restartstream >> c.parameters_t0[i] >> c.parameters_t1[i] >> c.dparameters_t0[i];
  }
  {
     auto & c = cFish->periodScheduler;
     restartstream >> c.t0 >> c.t1;
     for(int i=0;i<c.npoints;++i)
       restartstream >> c.parameters_t0[i] >> c.parameters_t1[i] >> c.dparameters_t0[i];
  }
  {
     auto & c = cFish->rlBendingScheduler;
     restartstream >> c.t0 >> c.t1;
     for(int i=0;i<c.npoints;++i)
       restartstream >> c.parameters_t0[i] >> c.parameters_t1[i] >> c.dparameters_t0[i];
  }
  {
     auto & c = cFish->torsionScheduler;
     restartstream >> c.t0 >> c.t1;
     for(int i=0;i<c.npoints;++i)
       restartstream >> c.parameters_t0[i] >> c.parameters_t1[i] >> c.dparameters_t0[i];
  }
  {
     auto & c = cFish->turnZScheduler;
     restartstream >> c.t0 >> c.t1;
     for(int i=0;i<c.npoints;++i)
       restartstream >> c.parameters_t0[i] >> c.parameters_t1[i] >> c.dparameters_t0[i];
  }
  restartstream.close();

  bool ret = true;
  double temp;
  int temp1;
  ret = ret && 1==fscanf(f, "origC_x: %le\n", &temp); origC[0] = temp;
  ret = ret && 1==fscanf(f, "origC_y: %le\n", &temp); origC[1] = temp;
  ret = ret && 1==fscanf(f, "origC_z: %le\n", &temp); origC[2] = temp;
  ret = ret && 1==fscanf(f, "errP: %le\n", &temp); cFish->errP = temp;
  ret = ret && 1==fscanf(f, "errD: %le\n", &temp); cFish->errD = temp;
  ret = ret && 1==fscanf(f, "curv_PID_fac             : %le\n", &temp); cFish->curv_PID_fac              = temp;
  ret = ret && 1==fscanf(f, "curv_PID_dif             : %le\n", &temp); cFish->curv_PID_dif              = temp;
  ret = ret && 1==fscanf(f, "avgDeltaY                : %le\n", &temp); cFish->avgDeltaY                 = temp;
  ret = ret && 1==fscanf(f, "avgDangle                : %le\n", &temp); cFish->avgDangle                 = temp;
  ret = ret && 1==fscanf(f, "avgAngVel                : %le\n", &temp); cFish->avgAngVel                 = temp;
  ret = ret && 1==fscanf(f, "lastTact                 : %le\n", &temp); cFish->lastTact                  = temp;
  ret = ret && 1==fscanf(f, "lastCurv                 : %le\n", &temp); cFish->lastCurv                  = temp;
  ret = ret && 1==fscanf(f, "oldrCurv                 : %le\n", &temp); cFish->oldrCurv                  = temp;
  ret = ret && 1==fscanf(f, "periodPIDval             : %le\n", &temp); cFish->periodPIDval              = temp;
  ret = ret && 1==fscanf(f, "periodPIDdif             : %le\n", &temp); cFish->periodPIDdif              = temp;
  ret = ret && 1==fscanf(f, "time0                    : %le\n", &temp); cFish->time0                     = temp;
  ret = ret && 1==fscanf(f, "timeshift                : %le\n", &temp); cFish->timeshift                 = temp;
  ret = ret && 1==fscanf(f, "lastTime                 : %le\n", &temp); cFish->lastTime                  = temp;
  ret = ret && 1==fscanf(f, "lastAvel                 : %le\n", &temp); cFish->lastAvel                  = temp;
  ret = ret && 1==fscanf(f, "Ttorsion_start           : %le\n", &temp); cFish->Ttorsion_start            = temp;
  ret = ret && 1==fscanf(f, "Tman_start               : %le\n", &temp); cFish->Tman_start                = temp;
  ret = ret && 1==fscanf(f, "Tman_finish              : %le\n", &temp); cFish->Tman_finish               = temp;
  ret = ret && 1==fscanf(f, "Lman;                    : %le\n", &temp); cFish->Lman                      = temp;
  ret = ret && 1==fscanf(f, "current_period           : %le\n", &temp); cFish->current_period            = temp;
  ret = ret && 1==fscanf(f, "next_period              : %le\n", &temp); cFish->next_period               = temp;
  ret = ret && 1==fscanf(f, "transition_start         : %le\n", &temp); cFish->transition_start          = temp;
  ret = ret && 1==fscanf(f, "transition_duration      : %le\n", &temp); cFish->transition_duration       = temp;
  ret = ret && 1==fscanf(f, "torsionValues[0]         : %le\n", &temp); cFish->torsionValues[0]          = temp;
  ret = ret && 1==fscanf(f, "torsionValues[1]         : %le\n", &temp); cFish->torsionValues[1]          = temp;
  ret = ret && 1==fscanf(f, "torsionValues[2]         : %le\n", &temp); cFish->torsionValues[2]          = temp;
  ret = ret && 1==fscanf(f, "torsionValues_previous[0]: %le\n", &temp); cFish->torsionValues_previous[0] = temp;
  ret = ret && 1==fscanf(f, "torsionValues_previous[1]: %le\n", &temp); cFish->torsionValues_previous[1] = temp;
  ret = ret && 1==fscanf(f, "torsionValues_previous[2]: %le\n", &temp); cFish->torsionValues_previous[2] = temp;
  ret = ret && 1==fscanf(f, "TperiodPID               : %d\n", &temp1); cFish->TperiodPID      = temp1;
  ret = ret && 1==fscanf(f, "control_torsion          : %d\n", &temp1); cFish->control_torsion = temp1;
  if( (not ret) ) {
    printf("Error reading restart file. Aborting...\n");
    fflush(0); abort();
  }
}

StefanFish::StefanFish(SimulationData & s, ArgumentParser&p) : Fish(s, p)
{
  const Real Tperiod = p("-T").asDouble(1.0);
  const Real phaseShift = p("-phi").asDouble(0.0);
  const Real ampFac = p("-amplitudeFactor").asDouble(1.0);
  bCorrectTrajectory  = p("-Correct" ).asBool(false);
  bCorrectTrajectoryZ = p("-CorrectZ").asBool(false);
  bCorrectPosition    = p("-CorrectPosition" ).asBool(false);
  bCorrectPositionZ   = p("-CorrectPositionZ").asBool(false);
  bCorrectRoll        = p("-CorrectRoll").asBool(false);

  myFish = new CurvatureDefinedFishData(length, Tperiod, phaseShift, sim.hmin, ampFac);

  std::string heightName = p("-heightProfile").asString("baseline");
  std::string  widthName = p( "-widthProfile").asString("baseline");
  MidlineShapes::computeWidthsHeights(heightName, widthName, length, 
		                      myFish->rS, myFish->height, myFish->width, myFish->Nm, sim.rank);
  origC[0] = position[0];
  origC[1] = position[1];
  origC[2] = position[2];

  if ( (bCorrectTrajectory || bCorrectTrajectoryZ || 
        bCorrectPosition   || bCorrectPositionZ || bCorrectRoll ) && std::fabs(quaternion[0]-1) > 1e-6)
  {
    std::cout << "PID controller only works for zero initial angles." << std::endl;
    MPI_Abort(sim.comm,1);
  }

  if(sim.rank==0) printf("nMidline=%d, length=%f, Tperiod=%f, phaseShift=%f\n",myFish->Nm, length, Tperiod, phaseShift);
}

void StefanFish::create()
{
  const Real q[4] = {quaternion[0],quaternion[1],quaternion[2],quaternion[3]}; 

  //const Real angle_roll  = atan2(2.0 * (q[3] * q[2] + q[0] * q[1]) ,   1.0 - 2.0 * (q[1] * q[1] + q[2] * q[2]));
  const Real angle_pitch = asin (2.0 * (q[2] * q[0] - q[3] * q[1]));
  const Real angle_yaw   = atan2(2.0 * (q[3] * q[0] + q[1] * q[2]) , - 1.0 + 2.0 * (q[0] * q[0] + q[1] * q[1]));

  // Control pos diffs
  const double xDiff = (position[0] - origC[0])/length;
  const double yDiff = (position[1] - origC[1])/length;
  const Real relU = (transVel[0] + sim.uinf[0]) / length;
  const Real relV = (transVel[1] + sim.uinf[1]) / length;

  auto * const cFish = dynamic_cast<CurvatureDefinedFishData*>( myFish );

  const Real Tperiod = cFish->Tperiod;
  const Real angDiff = angle_yaw;// - origAng;
  const Real lastAngVel = cFish->lastAvel;
  // compute ang vel at t - 1/2 dt such that we have a better derivative:
  const Real aVelMidP = (angVel[2] + lastAngVel)*Tperiod/2;
  const Real aVelDiff = (angVel[2] - lastAngVel)*Tperiod/sim.dt;
  cFish->lastAvel = angVel[2]; // store for next time
  const Real DT = sim.dt/Tperiod;
  // derivatives of following 2 exponential averages:
  const Real velDAavg = (angDiff-cFish->avgDangle)/Tperiod + DT * angVel[2];
  const Real velDYavg = (  yDiff-cFish->avgDeltaY)/Tperiod + DT * relV;
  const Real velAVavg = 10*((aVelMidP-cFish->avgAngVel)/Tperiod +DT*aVelDiff);

  // exponential averages
  cFish->avgDangle = (1.0-   DT) * cFish->avgDangle +    DT * angDiff;
  cFish->avgDeltaY = (1.0-   DT) * cFish->avgDeltaY +    DT *   yDiff;
  cFish->avgAngVel = (1.0-10*DT) * cFish->avgAngVel + 10*DT *aVelMidP; //faster average
  const Real avgDangle = cFish->avgDangle;
  const Real avgDeltaY = cFish->avgDeltaY;

  //control yaw angle and positin in xy plane - not tested very well!
  if (bCorrectPosition && sim.dt>0)
  {
    // integral (averaged) and proportional absolute DY and their derivative
    const Real absPy = std::fabs(yDiff), absIy = std::fabs(avgDeltaY);
    const Real velAbsPy =     yDiff>0 ? relV     : -relV;
    const Real velAbsIy = avgDeltaY>0 ? velDYavg : -velDYavg;

    const Real IangPdy = avgDangle * absPy;
    const Real PangIdy = angDiff   * absIy;
    const Real IangIdy = avgDangle * absIy;
    //time derivatives of above three terms:
    const Real velIangPdy = velAbsPy * avgDangle + absPy * velDAavg;
    const Real velPangIdy = velAbsIy * angDiff   + absIy * angVel[2];
    const Real velIangIdy = velAbsIy * avgDangle + absIy * velDAavg;
    //If angle is positive: positive curvature only if Dy<0 (must go up)
    //If angle is negative: negative curvature only if Dy>0 (must go down)
    //when appropriate, coefficients are set to zero.
    const Real coefIangPdy = avgDangle *     yDiff < 0 ? 1 : 0;
    const Real coefPangIdy = angDiff   * avgDeltaY < 0 ? 1 : 0;
    const Real coefIangIdy = avgDangle * avgDeltaY < 0 ? 1 : 0;
    const Real totalTerm = coefIangPdy *    IangPdy + coefPangIdy *    PangIdy + coefIangIdy *    IangIdy;
    const Real totalDiff = coefIangPdy * velIangPdy + coefPangIdy * velPangIdy + coefIangIdy * velIangIdy;
    cFish->correctTrajectory(totalTerm, totalDiff);
    const Real periodFac = 1.0 - xDiff;
    const Real periodVel =     - relU;
    cFish->correctTailPeriod(periodFac, periodVel, sim.time, sim.dt);
  }
  //control yaw angle
  else if (bCorrectTrajectory && sim.dt>0)
  {
    const Real avgAngVel = cFish->avgAngVel, absAngVel = std::fabs(avgAngVel);
    const Real absAvelDiff = avgAngVel>0? velAVavg : -velAVavg;
    const Real coefInst = angDiff*avgAngVel>0 ? 0.01 : 1, coefAvg = 0.1;
    const Real termInst = angDiff*absAngVel;
    const Real diffInst = angDiff*absAvelDiff + angVel[2]*absAngVel;
    Real totalTerm = coefInst*termInst + coefAvg*avgDangle;
    Real totalDiff = coefInst*diffInst + coefAvg*velDAavg;
    if (totalTerm >  0.2) {totalTerm =  0.2; totalDiff = 0;}
    if (totalTerm < -0.2) {totalTerm = -0.2; totalDiff = 0;}
    cFish->correctTrajectory(totalTerm, totalDiff);
  }

  //const Real periodFac = 1.0 - xDiff;
  //const Real periodVel =     - relU;
  //cFish->correctTailPeriod(periodFac, periodVel, sim.time, sim.dt);
  //if (sim.rank == 0 && sim.step % 1 == 0)
  //{
  //  char buf[500];
  //  sprintf(buf, "angles_%d.dat", obstacleID);
  //  FILE * f = fopen(buf, "a");
  //  fprintf(f, "%g %g %g %g \n",sim.time,angle_roll*180/M_PI,angle_pitch*180/M_PI,angle_yaw*180/M_PI);
  //  fclose(f);
  //}

  //control pitching
  if (!bCorrectPositionZ && bCorrectTrajectoryZ && sim.dt > 0)
  {
    const Real dqdt[4] = { 0.5*( - angVel[0]*q[1] - angVel[1]*q[2] - angVel[2]*q[3] ),
                           0.5*( + angVel[0]*q[0] + angVel[1]*q[3] - angVel[2]*q[2] ),
                           0.5*( - angVel[0]*q[3] + angVel[1]*q[0] + angVel[2]*q[1] ),
                           0.5*( + angVel[0]*q[2] - angVel[1]*q[1] + angVel[2]*q[0] )};

    const Real arg_aux = 2.0 * (q[2] * q[0] - q[3] * q[1]);

    const Real dpitch_dt = 2.0 / ( sqrt(1.0 - arg_aux*arg_aux) + 1e-21 ) * (q[2]*dqdt[0]+dqdt[2]*q[0]-q[3]*dqdt[1]-dqdt[3]*q[1]);

    const Real rel = min(1.,100*sim.dt/Tperiod);

    cFish->errP = - angle_pitch;
    cFish->errD = (1-rel) * cFish->errD + rel * (-dpitch_dt);

    const Real u = 4.0*cFish->errP + 5.0*cFish->errD;
    //if (std::fabs(u) > 5.) u = u > 0 ? 5. : -5.;
    cFish->action_torsion_pitching_radius(sim.time, sim.time, -u/length);
  }
  //control position in Z plane and pitching with PD controller - not tested very well!
  else if (bCorrectPositionZ) 
  {
    const Real dqdt[4] = { 0.5*( - angVel[0]*q[1] - angVel[1]*q[2] - angVel[2]*q[3] ),
                           0.5*( + angVel[0]*q[0] + angVel[1]*q[3] - angVel[2]*q[2] ),
                           0.5*( - angVel[0]*q[3] + angVel[1]*q[0] + angVel[2]*q[1] ),
                           0.5*( + angVel[0]*q[2] - angVel[1]*q[1] + angVel[2]*q[0] )};

    const Real arg_aux = 2.0 * (q[2] * q[0] - q[3] * q[1]);

    const Real dpitch_dt = 2.0 / ( sqrt(1.0 - arg_aux*arg_aux) + 1e-21 ) * (q[2]*dqdt[0]+dqdt[2]*q[0]-q[3]*dqdt[1]-dqdt[3]*q[1]);

    const Real rel = min(1.,100*sim.dt/Tperiod);

    cFish->errP = - angle_pitch + 0.1*(origC[2] - position[2])/length;
    cFish->errD = (1-rel) * cFish->errD + rel * (-dpitch_dt - 0.1*transVel[2]/length);

    const Real u = 4.0*cFish->errP + 5.0*cFish->errD;
    cFish->action_torsion_pitching_radius(sim.time, sim.time, -u/length);
  }

  Fish::create();
}

void StefanFish::computeVelocities()
{
  Obstacle::computeVelocities();

  if (bCorrectRoll) //kill angVel component the causes rolling
  {
    const Real q[4] = {quaternion[0],quaternion[1],quaternion[2],quaternion[3]}; 
    const Real angle_roll  = atan2(2.0 * (q[3] * q[2] + q[0] * q[1]) ,   1.0 - 2.0 * (q[1] * q[1] + q[2] * q[2]));

    const Real R[3] = { 1.0 - 2*(q[2]*q[2]+q[3]*q[3]), 
                                2*(q[1]*q[2]+q[3]*q[0]), 
                                2*(q[1]*q[3]-q[2]*q[0])}; 
    //const Real R[3] = {       2*(q[1]*q[2]-q[3]*q[0]), 
    //                      1.0 - 2*(q[1]*q[1]+q[3]*q[3]),
    //                            2*(q[2]*q[3]+q[1]*q[0])}; 
    //const Real R[3] = {       2*(q[1]*q[3]+q[2]*q[0]), 
    //                            2*(q[2]*q[3]-q[1]*q[0]), 
    //                      1.0 - 2*(q[1]*q[1]+q[2]*q[2])};
    const Real unit_vector[3] = {R[0],R[1],R[2]};

    //kill angVel in unit_vector direction
    const Real mag = angVel[0]*unit_vector[0] + angVel[1]*unit_vector[1] + angVel[2]*unit_vector[2];
    const Real angVel_projection[3] = {mag*unit_vector[0],mag*unit_vector[1],mag*unit_vector[2]};
    angVel[0] = angVel[0] - angVel_projection[0] - angle_roll/sim.dt*unit_vector[0];
    angVel[1] = angVel[1] - angVel_projection[1] - angle_roll/sim.dt*unit_vector[1];
    angVel[2] = angVel[2] - angVel_projection[2] - angle_roll/sim.dt*unit_vector[2];
  }
}

//////////////////////////////////
//Reinforcement Learning functions
//////////////////////////////////

void StefanFish::act(const Real t_rlAction, const std::vector<Real>& a) const
{
  auto * const cFish = dynamic_cast<CurvatureDefinedFishData*>( myFish );
  if( cFish == nullptr ) { printf("Someone touched my fish\n"); abort(); }
  std::vector <Real> actions = a;
  if (actions.size() == 0)
  {
     std::cerr << "No actions given to CurvatureDefinedFishData::execute\n";
     MPI_Abort(sim.comm,1);
  }
  if (bForcedInSimFrame[2] == true && a.size() > 1) actions[1] = 0; //no pitching
  cFish->execute(sim.time, t_rlAction, actions);
}

Real StefanFish::getLearnTPeriod() const
{
  auto * const cFish = dynamic_cast<CurvatureDefinedFishData*>( myFish );
  assert( cFish != nullptr);
  return cFish->next_period;
}

Real StefanFish::getPhase(const Real t) const
{
  auto * const cFish = dynamic_cast<CurvatureDefinedFishData*>( myFish );
  const Real T0 = cFish->time0;
  const Real Ts = cFish->timeshift;
  const Real Tp = cFish->periodPIDval;
  const Real arg  = 2*M_PI*((t-T0)/Tp +Ts) + M_PI*cFish->phaseShift;
  const Real phase = std::fmod(arg, 2*M_PI);
  return (phase<0) ? 2*M_PI + phase : phase;
}

std::vector<Real> StefanFish::state() const
{
  auto * const cFish = dynamic_cast<CurvatureDefinedFishData*>( myFish );
  assert( cFish != nullptr);
  const Real Tperiod = cFish->Tperiod;
  std::vector<Real> S(25);
  S[0 ] = position[0];
  S[1 ] = position[1];
  S[2 ] = position[2];

  S[3 ] = quaternion[0];
  S[4 ] = quaternion[1];
  S[5 ] = quaternion[2];
  S[6 ] = quaternion[3];

  S[7 ] = getPhase(sim.time);

  S[8 ] = transVel[0] * Tperiod / length;
  S[9 ] = transVel[1] * Tperiod / length;
  S[10] = transVel[2] * Tperiod / length;

  S[11] = angVel[0] * Tperiod;
  S[12] = angVel[1] * Tperiod;
  S[13] = angVel[2] * Tperiod;

  S[14] = cFish->lastCurv;
  S[15] = cFish->oldrCurv;

  //sensor locations
  const std::array<Real,3> locFront = {cFish->sensorLocation[0*3+0],cFish->sensorLocation[0*3+1],cFish->sensorLocation[0*3+2]};
  const std::array<Real,3> locUpper = {cFish->sensorLocation[1*3+0],cFish->sensorLocation[1*3+1],cFish->sensorLocation[1*3+2]};
  const std::array<Real,3> locLower = {cFish->sensorLocation[2*3+0],cFish->sensorLocation[2*3+1],cFish->sensorLocation[2*3+2]};
  //compute shear stress force (x,y,z) components
  std::array<Real,3> shearFront = getShear( locFront );
  std::array<Real,3> shearUpper = getShear( locLower );
  std::array<Real,3> shearLower = getShear( locUpper );
  S[16] = shearFront[0] * Tperiod / length;
  S[17] = shearFront[1] * Tperiod / length;
  S[18] = shearFront[2] * Tperiod / length;
  S[19] = shearUpper[0] * Tperiod / length;
  S[20] = shearUpper[1] * Tperiod / length;
  S[21] = shearUpper[2] * Tperiod / length;
  S[22] = shearLower[0] * Tperiod / length;
  S[23] = shearLower[1] * Tperiod / length;
  S[24] = shearLower[2] * Tperiod / length;
  return S;
}

ssize_t StefanFish::holdingBlockID(const std::array<Real,3> pos) const
{
  const std::vector<cubism::BlockInfo>& velInfo = sim.velInfo();
  for(size_t i=0; i<velInfo.size(); ++i)
  {
    // compute lower left and top right corners of block (+- 0.5 h because pos returns cell centers)
    std::array<Real,3> MIN = velInfo[i].pos<Real>(0                  , 0                 ,  0                  );
    std::array<Real,3> MAX = velInfo[i].pos<Real>(ScalarBlock::sizeX-1, ScalarBlock::sizeY-1, ScalarBlock::sizeZ-1);
    MIN[0] -= 0.5 * velInfo[i].h;
    MIN[1] -= 0.5 * velInfo[i].h;
    MIN[2] -= 0.5 * velInfo[i].h;
    MAX[0] += 0.5 * velInfo[i].h;
    MAX[1] += 0.5 * velInfo[i].h;
    MAX[2] += 0.5 * velInfo[i].h;

    // check whether point is inside block
    if( pos[0] >= MIN[0] && pos[1] >= MIN[1] && pos[2] >= MIN[2] && pos[0] <= MAX[0] && pos[1] <= MAX[1] && pos[2] <= MAX[2] )
    {
      return i;
    }
  }
  return -1; // rank does not contain point
};

// returns shear at given surface location
std::array<Real, 3> StefanFish::getShear(const std::array<Real,3> pSurf) const
{
  const std::vector<cubism::BlockInfo>& velInfo = sim.velInfo(); 

  Real myF[3] = {0,0,0};
  
  // Get blockId of block that contains point pSurf.
  ssize_t blockIdSurf = holdingBlockID(pSurf);
  if( blockIdSurf >= 0 )
  {
    const auto & skinBinfo = velInfo[blockIdSurf];

    // check whether obstacle block exists
    if(obstacleBlocks[blockIdSurf] != nullptr )
    {
      Real dmin = 1e10;
      ObstacleBlock * const O = obstacleBlocks[blockIdSurf];
      for(int k = 0; k < O->nPoints; ++k)
      {
        const int ix = O->surface[k]->ix;
        const int iy = O->surface[k]->iy;
        const int iz = O->surface[k]->iz;
        const std::array<Real,3> p = skinBinfo.pos<Real>(ix, iy, iz);
        const Real d = (p[0]-pSurf[0])*(p[0]-pSurf[0])+(p[1]-pSurf[1])*(p[1]-pSurf[1])+(p[2]-pSurf[2])*(p[2]-pSurf[2]);
        if (d < dmin)
        {
          dmin = d;
          myF[0] = O->fxV[k];
          myF[1] = O->fyV[k];
          myF[2] = O->fzV[k];
        }
      }
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, myF, 3, MPI_Real, MPI_SUM, sim.comm);

  return std::array<Real, 3>{{myF[0],myF[1],myF[2]}};// return shear
};

CubismUP_3D_NAMESPACE_END
