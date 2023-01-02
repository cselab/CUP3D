//
//  Cubism3D
//  Copyright (c) 2022 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//

#include "StefanFish.h"

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

  const std::array<Real,6> curvaturePoints = { 0.0,0.15*length,0.4*length,0.65*length,0.9*length,length};
  const std::array<Real,7>      bendPoints = { -0.5,-0.25,0.0,0.25,0.5,0.75,1.0};
  const std::array<Real,6> curvatureValues = { 0.82014/length, 1.46515/length, 2.57136/length, 3.75425/length, 5.09147/length, 5.70449/length};

  #if 1 // ramp-up over Tperiod
  const std::array<Real,6> curvatureZeros = std::array<Real, 6>();
  curvatureScheduler.transition(0,0,Tperiod,curvatureZeros,curvatureValues);
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
    const Real  curv = std::sin(arg)     + rB[i] +  beta;
    const Real dcurv = std::cos(arg)*darg+ vB[i] + dbeta;
    rK[i] = alpha*amplitudeFactor* rC[i]*curv;
    vK[i] = alpha*amplitudeFactor*(vC[i]*curv + rC[i]*dcurv)+ dalpha*amplitudeFactor*rC[i]*curv;
    rT[i] = 0;
    vT[i] = 0;
    assert(!std::isnan(rK[i]));
    assert(!std::isinf(rK[i]));
    assert(!std::isnan(vK[i]));
    assert(!std::isinf(vK[i]));
  }

  if (control_torsion)
  {
    const std::array<Real ,3> torsionPoints = { 0.0, 0.5*length, length };
    torsionScheduler.transition (t,Ttorsion_start,Ttorsion_start+0.5*Tperiod,torsionValues_previous,torsionValues);
    torsionScheduler.gimmeValues(t,torsionPoints,Nm,rS,rT,vT);
  }
  Frenet3D::solve(Nm, rS, rK, vK, rT, vT, rX, rY, rZ, vX, vY, vZ, norX, norY, norZ, vNorX, vNorY, vNorZ, binX, binY, binZ, vBinX, vBinY, vBinZ);

  performPitchingMotion(t);
}

void CurvatureDefinedFishData::performPitchingMotion(const Real t)
{
  Real R,Rdot;

  if (std::fabs(gamma) > 1e-10) {R = 1.0/gamma; Rdot = - 1.0/gamma/gamma * dgamma;}
  else {R= gamma >= 0 ? 1e10:-1e10; Rdot=0.0;}

  const Real x0N    = rX[Nm-1];
  const Real y0N    = rY[Nm-1];
  const Real x0Ndot = vX[Nm-1];
  const Real y0Ndot = vY[Nm-1];
  const Real phi    = atan2(y0N,x0N);
  const Real phidot = 1.0/(1.0+pow(y0N/x0N,2))*(y0Ndot/x0N-y0N*x0Ndot/x0N/x0N);
  const Real M      = pow(x0N*x0N+y0N*y0N,0.5);
  const Real Mdot   = (x0N*x0Ndot+y0N*y0Ndot)/M;
  const Real cosphi = cos(phi);
  const Real sinphi = sin(phi);
  #pragma omp parallel for
  for(int i=0; i<Nm; i++)
  {
    const double x0    = rX[i];
    const double y0    = rY[i];
    const double x0dot = vX[i];
    const double y0dot = vY[i];
    const double x1    = cosphi*x0-sinphi*y0;
    const double y1    = sinphi*x0+cosphi*y0;
    const double x1dot = cosphi*x0dot-sinphi*y0dot+(-sinphi*x0-cosphi*y0)*phidot;
    const double y1dot = sinphi*x0dot+cosphi*y0dot+( cosphi*x0-sinphi*y0)*phidot;
    const double theta = (M-x1)/R;
    const double costheta = cos(theta);
    const double sintheta = sin(theta);
    const double x2       = M    - R * sintheta;
    const double y2       = y1;
    const double z2       = R - R * costheta;
    const double thetadot = (Mdot-x1dot)/R-(M-x1)/R/R*Rdot;
    const double x2dot = Mdot - Rdot*sintheta-R*costheta*thetadot;
    const double y2dot = y1dot;
    const double z2dot = Rdot-Rdot*costheta+R*sintheta*thetadot;
    rX[i] = x2;
    rY[i] = y2;
    rZ[i] = z2;
    vX[i] = x2dot;
    vY[i] = y2dot;
    vZ[i] = z2dot;
  }

  recomputeNormalVectors();
}

void CurvatureDefinedFishData::recomputeNormalVectors()
{
  //compute normal and binormal vectors for a given midline
  #pragma omp parallel for
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
        savestream << c.parameters_t0[i]  << "\t" << c.parameters_t1[i]  << "\t" << c.dparameters_t0[i] << std::endl;
    }
    {
      const auto & c = cFish->periodScheduler;
      savestream << c.t0 << "\t" << c.t1 << std::endl;
      for(int i=0;i<c.npoints;++i)
        savestream << c.parameters_t0[i]  << "\t" << c.parameters_t1[i]  << "\t" << c.dparameters_t0[i] << std::endl;
    }
    {
      const auto & c = cFish->rlBendingScheduler;
      savestream << c.t0 << "\t" << c.t1 << std::endl;
      for(int i=0;i<c.npoints;++i)
          savestream << c.parameters_t0[i]  << "\t" << c.parameters_t1[i]  << "\t" << c.dparameters_t0[i] << std::endl;
    }
    {
      const auto & c = cFish->torsionScheduler;
      savestream << c.t0 << "\t" << c.t1 << std::endl;
      for(int i=0;i<c.npoints;++i)
         savestream << c.parameters_t0[i]  << "\t" << c.parameters_t1[i]  << "\t" << c.dparameters_t0[i] << std::endl;
    }
    savestream.close();
  }

  //Save these numbers for PID controller and other stuff. Maybe not all of them are needed
  //but we don't care, it's only a few numbers.
  fprintf(f, "origC_x: %20.20e\n", (double)origC[0]);
  fprintf(f, "origC_y: %20.20e\n", (double)origC[1]);
  fprintf(f, "origC_z: %20.20e\n", (double)origC[2]);
  fprintf(f,"lastTact                 : %20.20e\n",(double)cFish->lastTact                 );
  fprintf(f,"lastCurv                 : %20.20e\n",(double)cFish->lastCurv                 );
  fprintf(f,"oldrCurv                 : %20.20e\n",(double)cFish->oldrCurv                 );
  fprintf(f,"periodPIDval             : %20.20e\n",(double)cFish->periodPIDval             );
  fprintf(f,"periodPIDdif             : %20.20e\n",(double)cFish->periodPIDdif             );
  fprintf(f,"time0                    : %20.20e\n",(double)cFish->time0                    );
  fprintf(f,"timeshift                : %20.20e\n",(double)cFish->timeshift                );
  fprintf(f,"Ttorsion_start           : %20.20e\n",(double)cFish->Ttorsion_start           );
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
  fprintf(f,"alpha                    : %20.20e\n",(double)cFish->alpha                    );
  fprintf(f,"dalpha                   : %20.20e\n",(double)cFish->dalpha                   );
  fprintf(f,"beta                     : %20.20e\n",(double)cFish->beta                     );
  fprintf(f,"dbeta                    : %20.20e\n",(double)cFish->dbeta                    );
  fprintf(f,"gamma                    : %20.20e\n",(double)cFish->gamma                    );
  fprintf(f,"dgamma                   : %20.20e\n",(double)cFish->dgamma                   );
  //TODO:save vector of integral terms for error in PID
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
  restartstream.close();
 
  bool ret = true;
  double temp;
  int temp1;
  ret = ret && 1==fscanf(f, "origC_x: %le\n", &temp); origC[0] = temp;
  ret = ret && 1==fscanf(f, "origC_y: %le\n", &temp); origC[1] = temp;
  ret = ret && 1==fscanf(f, "origC_z: %le\n", &temp); origC[2] = temp;
  ret = ret && 1==fscanf(f, "lastTact                 : %le\n", &temp); cFish->lastTact                  = temp;
  ret = ret && 1==fscanf(f, "lastCurv                 : %le\n", &temp); cFish->lastCurv                  = temp;
  ret = ret && 1==fscanf(f, "oldrCurv                 : %le\n", &temp); cFish->oldrCurv                  = temp;
  ret = ret && 1==fscanf(f, "periodPIDval             : %le\n", &temp); cFish->periodPIDval              = temp;
  ret = ret && 1==fscanf(f, "periodPIDdif             : %le\n", &temp); cFish->periodPIDdif              = temp;
  ret = ret && 1==fscanf(f, "time0                    : %le\n", &temp); cFish->time0                     = temp;
  ret = ret && 1==fscanf(f, "timeshift                : %le\n", &temp); cFish->timeshift                 = temp;
  ret = ret && 1==fscanf(f, "Ttorsion_start           : %le\n", &temp); cFish->Ttorsion_start            = temp;
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
  ret = ret && 1==fscanf(f, "alpha                    : %le\n", &temp); cFish->alpha                     = temp;
  ret = ret && 1==fscanf(f, "dalpha                   : %le\n", &temp); cFish->dalpha                    = temp;
  ret = ret && 1==fscanf(f, "beta                     : %le\n", &temp); cFish->beta                      = temp;
  ret = ret && 1==fscanf(f, "dbeta                    : %le\n", &temp); cFish->dbeta                     = temp;
  ret = ret && 1==fscanf(f, "gamma                    : %le\n", &temp); cFish->gamma                     = temp;
  ret = ret && 1==fscanf(f, "dgamma                   : %le\n", &temp); cFish->dgamma                    = temp;

  if( (not ret) )
   {
    printf("Error reading restart file. Aborting...\n");
    fflush(0); abort();
  }
}

StefanFish::StefanFish(SimulationData & s, ArgumentParser&p) : Fish(s, p)
{
  const Real Tperiod     = p("-T"               ).asDouble(1.0);
  const Real phaseShift  = p("-phi"             ).asDouble(0.0);
  const Real ampFac      = p("-amplitudeFactor" ).asDouble(1.0);
  bCorrectPosition       = p("-CorrectPosition" ).asBool(false);
  bCorrectPositionZ      = p("-CorrectPositionZ").asBool(false);
  bCorrectRoll           = p("-CorrectRoll"     ).asBool(false);
  std::string heightName = p("-heightProfile"   ).asString("baseline");
  std::string  widthName = p( "-widthProfile"   ).asString("baseline");

  if ( (bCorrectPosition   || bCorrectPositionZ || bCorrectRoll ) && std::fabs(quaternion[0]-1) > 1e-6)
  {
    std::cout << "PID controller only works for zero initial angles." << std::endl;
    MPI_Abort(sim.comm,1);
  }

  myFish = new CurvatureDefinedFishData(length, Tperiod, phaseShift, sim.hmin, ampFac);

  MidlineShapes::computeWidthsHeights(heightName, widthName, length, myFish->rS, myFish->height, myFish->width, myFish->Nm, sim.rank);
  origC[0] = position[0];
  origC[1] = position[1];
  origC[2] = position[2];

  if(sim.rank==0) printf("nMidline=%d, length=%f, Tperiod=%f, phaseShift=%f\n",myFish->Nm, length, Tperiod, phaseShift);

  wyp = p("-wyp").asDouble(1.0);
  wzp = p("-wzp").asDouble(1.0);
}

static void clip_quantities(const Real fmax, const Real dfmax, const Real dt, const bool zero, const Real fcandidate, const Real dfcandidate, Real & f, Real & df)
{
    if (zero)
    {
      f = 0;
      df = 0;
    }
    else if (std::fabs(dfcandidate) > dfmax)
    {
      df = dfcandidate > 0 ? +dfmax : -dfmax;
      f  = f + dt * df;
    }
    else if (std::fabs(fcandidate) < fmax)
    {
      f  =  fcandidate;
      df = dfcandidate;
    }
    else
    {
      f = fcandidate > 0 ? fmax : -fmax;
      df = 0;
    }
}

void StefanFish::create()
{
  const Real q [4] = {quaternion[0],quaternion[1],quaternion[2],quaternion[3]};
  const Real dq[4] = { 0.5*( - angVel[0]*q[1] - angVel[1]*q[2] - angVel[2]*q[3] ),
                       0.5*( + angVel[0]*q[0] + angVel[1]*q[3] - angVel[2]*q[2] ),
                       0.5*( - angVel[0]*q[3] + angVel[1]*q[0] + angVel[2]*q[1] ),
                       0.5*( + angVel[0]*q[2] - angVel[1]*q[1] + angVel[2]*q[0] )};
  auto * const cFish = dynamic_cast<CurvatureDefinedFishData*>( myFish );
  const Real angle_roll  = atan2(2.0 * (q[3] * q[2] + q[0] * q[1]) ,   1.0 - 2.0 * (q[1] * q[1] + q[2] * q[2]));
  const bool roll_is_small = std::fabs(angle_roll) < 20* M_PI/180.;

  if (bCorrectPosition)
  {
    //1.control position in x
    cFish->alpha  = 1.0 + (position[0]               - origC[0])/length;
    cFish->dalpha =       (transVel[0] + sim.uinf[0]           )/length;
    if (roll_is_small == false)
    {
      cFish-> alpha = 1.0;
      cFish->dalpha = 0.0;
    }

    //2.control position in y and yaw angle
    const Real y        = absPos[1];
    const Real ytgt     =  origC[1];
    const Real dy       = (ytgt-y          )/length;
    const Real dydt     = (    -transVel[1])/length;
    const Real signY    = dy > 0 ? 1 : -1;
    const Real  yaw_tgt = 0;
    const Real dyaw_tgt = 0;
    const Real    nom   =         2.0 * (q[3] * q[0] + q[1] * q[2]);
    const Real  denom   = - 1.0 + 2.0 * (q[0] * q[0] + q[1] * q[1]);
    const Real    yaw   = atan2( nom , denom);
    const Real   dphi   =  yaw- yaw_tgt;
    const Real   dnom   = 2.0 * (dq[3] * q[0] + dq[0] * q[3]+q[1] * dq[2] + q[2] * dq[1]);
    const Real ddenom   = 2.0 * (2.0*q[0] * dq[0] + 2.0*q[1] * dq[1]);
    const Real    arg   = nom/denom;
    const Real   darg   = (dnom*denom-nom*ddenom)/denom/denom;
    const Real   dyaw   = 1.0/(1.0+arg*arg)*darg; 
    const Real dphidt   = dyaw-dyaw_tgt;
    const Real b        = wyp * signY * dy * dphi;
    const Real dbdt     = wyp * signY *(dydt*dphi + dy*dphidt);
    clip_quantities(0.5,0.1,sim.dt,!roll_is_small,b,dbdt,cFish->beta,cFish->dbeta);
  }
  if (bCorrectPositionZ)
  {
    //compute pitch and d(pitch)/dt
    Real pitch,dpitch;
    {
      //pitch = asin (2.0 * (q[2] * q[0] - q[3] * q[1]));
      //const Real arg_aux = 2.0 * (q[2] * q[0] - q[3] * q[1]);
      //dpitch = 2.0 / ( sqrt(1.0 - arg_aux*arg_aux) + 1e-21 ) * (q[2]*dq[0]+dq    [2]*q[0]-q[3]*dq[1]-dq[3]*q[1]);

      const int  Nm = cFish->Nm;

      const Real Rmatrix3D[3] = {2*(q[1]*q[3]-q[2]*q[0]), 
                                 2*(q[2]*q[3]+q[1]*q[0]),
                               1-2*(q[1]*q[1]+q[2]*q[2])};
      const Real d1 = cFish->rX[0]-cFish->rX[Nm/2];
      const Real d2 = cFish->rY[0]-cFish->rY[Nm/2];
      const Real d3 = cFish->rZ[0]-cFish->rZ[Nm/2];
      const Real dn = pow(d1*d1+d2*d2+d3*d3,0.5)+1e-21;
      const Real vx = d1/dn;
      const Real vy = d2/dn;
      const Real vz = d3/dn;
      const Real xx2 = Rmatrix3D[0]*vx + Rmatrix3D[1]*vy + Rmatrix3D[2]*vz;
      pitch = asin(xx2);

      const Real dR[3] = {2*(dq[1]*q[3]-dq[2]*q[0] + q[1]*dq[3]-q[2]*dq[0]), 
                          2*(dq[2]*q[3]+dq[1]*q[0] + q[2]*dq[3]+q[1]*dq[0]), 
                         -2*(dq[1]*q[1]+dq[2]*q[2] + q[1]*dq[1]+q[2]*dq[2])};
      const Real dd1 = cFish->vX[0]-cFish->vX[Nm/2];
      const Real dd2 = cFish->vY[0]-cFish->vY[Nm/2];
      const Real dd3 = cFish->vZ[0]-cFish->vZ[Nm/2];
      const Real ddn = pow(d1*d1+d2*d2+d3*d3,-0.5)*(d1*dd1+d2*dd2+d3*dd3);
      const Real dvx = (dd1*dn-d1*ddn)/(dn*dn);
      const Real dvy = (dd2*dn-d2*ddn)/(dn*dn);
      const Real dvz = (dd3*dn-d3*ddn)/(dn*dn);
      const Real dxx2 = dR[0]*vx+dR[1]*vy+dR[2]*vz+ Rmatrix3D[0]*dvx + Rmatrix3D[1]*dvy + Rmatrix3D[2]*dvz;
      dpitch = 1.0 / ( sqrt(1.0 - xx2*xx2) + 1e-21) * dxx2;
    }

    const Real z          = absPos[2];
    const Real ztgt       = origC[2];
    const Real  pitch_tgt = 0;
    const Real dpitch_tgt = 0;
    const Real dz         = (ztgt-z          )/length;
    const Real dzdt       = (    -transVel[2])/length;
    const Real dphi       =  pitch- pitch_tgt;
    const Real dphidt     = dpitch-dpitch_tgt;
    const Real signZ      = dz > 0 ? 1 : -1;
    const Real g          = -wzp * signZ * dz * dphi;
    const Real dgdt       = -wzp * signZ * (dzdt * dphi + dz * dphidt);
    const Real gmax       = 1.0/length;
    const Real dRdtmax    = length/cFish->Tperiod;
    const Real dgdtmax    = std::fabs(gmax*gmax*dRdtmax);
    clip_quantities(gmax,dgdtmax,sim.dt,false,g,dgdt,cFish->gamma,cFish->dgamma);
  }

  #if 0
  if (sim.rank == 0)
  {
      char buf[500];
      sprintf(buf, "gamma%d.txt",obstacleID);
      FILE * f = fopen(buf, "a");
      fprintf(f, "%g %g %g %g %g \n",sim.time,cFish->beta,cFish->dbeta,cFish->gamma,cFish->dgamma);
      fclose(f);
  }
  #endif

  Fish::create();
}

void StefanFish::computeVelocities()
{
  Obstacle::computeVelocities();

  //Compute angular velocity component on the rolling axis of the fish and set it to 0.
  //Then, impose rolling angular velocity that will make the rolling angle go to zero after 
  //0.5Tperiod time has passed.
  //Important: this assumes an initial orientation (q = (1,0,0,0)) along the x-axis for the fish
  //           where the head is facing the (-1,0,0) direction 
  //           (this is the default initial orientation for fish). 
  if (bCorrectRoll)
  {
    auto * const cFish = dynamic_cast<CurvatureDefinedFishData*>( myFish );
    const Real T = cFish->Tperiod;
    const Real q[4] = {quaternion[0],quaternion[1],quaternion[2],quaternion[3]}; 
    const Real angle_roll  = atan2(2.0 * (q[3] * q[2] + q[0] * q[1]) ,   1.0 - 2.0 * (q[1] * q[1] + q[2] * q[2]));

    //const Real dq[4] = { 0.5*( - angVel[0]*q[1] - angVel[1]*q[2] - angVel[2]*q[3] ),
    //                     0.5*( + angVel[0]*q[0] + angVel[1]*q[3] - angVel[2]*q[2] ),
    //                     0.5*( - angVel[0]*q[3] + angVel[1]*q[0] + angVel[2]*q[1] ),
    //                     0.5*( + angVel[0]*q[2] - angVel[1]*q[1] + angVel[2]*q[0] )};
    //const Real  nom = 2.0 * (q[3] * q[2] + q[0] * q[1]);
    //const Real dnom = 2.0 * (dq[3] * q[2] + dq[0] * q[1]+q[3] * dq[2] + q[0] * dq[1]);
    //const Real denom = 1.0 - 2.0 * (q[1] * q[1] + q[2] * q[2]);
    //const Real ddenom = - 2.0 * (2.0*q[1] * dq[1] + 2.0*q[2] * dq[2]);
    //const Real  arg = nom/denom;
    //const Real darg = (dnom*denom-nom*ddenom)/denom/denom;
    //const Real da = 1.0/(1.0+arg*arg)*darg; 

    //1st column of rotation matrix
    const Real R0[3] = { 1.0 - 2*(q[2]*q[2]+q[3]*q[3]), 
                               2*(q[1]*q[2]+q[3]*q[0]), 
                               2*(q[1]*q[3]-q[2]*q[0])}; 
    //const Real R1[3] = {       2*(q[1]*q[2]-q[3]*q[0]), 
    //                     1.0 - 2*(q[1]*q[1]+q[3]*q[3]),
    //                           2*(q[2]*q[3]+q[1]*q[0])}; 
    //const Real R2[3] = {       2*(q[1]*q[3]+q[2]*q[0]), 
    //                           2*(q[2]*q[3]-q[1]*q[0]), 
    //                     1.0 - 2*(q[1]*q[1]+q[2]*q[2])};

    //current orientation is R * (1,0,0), where (1,0,0) is the initial (assumed) orientation
    //with which we only need the first column of R.

    const int sign_roll       = angle_roll > 0 ? 1:-1;
    const Real unit_vector[3] = {sign_roll*R0[0],sign_roll*R0[1],sign_roll*R0[2]};

    //const Real angVel_roll = unit_vector[0]*angVel[0]+unit_vector[1]*angVel[1]+unit_vector[2]*angVel[2];
    //if (angVel_roll < 0) return;

    const Real correction_magnitude = std::fabs(angle_roll)/(0.1*T);
    angVel[0] -= correction_magnitude*unit_vector[0];
    angVel[1] -= correction_magnitude*unit_vector[1];
    angVel[2] -= correction_magnitude*unit_vector[2];
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
