//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#include "Naca.h"

#include <cmath>

CubismUP_3D_NAMESPACE_BEGIN
using namespace cubism;

class NacaMidlineData : public FishMidlineData
{
 Real * const rK;
 Real * const vK;
 Real * const rC;
 Real * const vC;
 public:
  NacaMidlineData(const Real L, const Real _h, Real zExtent,
    Real t_ratio, Real HoverL=1) : FishMidlineData(L, 1, 0, _h),
    rK(_alloc(Nm)),vK(_alloc(Nm)), rC(_alloc(Nm)),vC(_alloc(Nm))
  {
    for(int i=0;i<Nm;++i) height[i] = length*HoverL/2;
    MidlineShapes::naca_width(t_ratio, length, rS, width, Nm);

    computeMidline(0.0, 0.0);
  }

  void computeMidline(const Real time, const Real dt) override
  {
    #if 1
      rX[0] = rY[0] = rZ[0] = 0.0;
      vX[0] = vY[0] = vZ[0] = 0.0;
      norX[0] = 0.0; norY[0] = 1.0; norZ[0] = 0.0;
      binX[0] = 0.0; binY[0] = 0.0; binZ[0] = 1.0;
      vNorX[0] = vNorY[0] = vNorZ[0] = 0.0;
      vBinX[0] = vBinY[0] = vBinZ[0] = 0.0;
      for(int i=1; i<Nm; ++i)
      {
        rY[i] = rZ[i] = 0.0;
        vX[i] = vY[i] = vZ[i] = 0.0;
        rX[i] = rX[i-1] + std::fabs(rS[i]-rS[i-1]);
        norX[i] = 0.0; norY[i] = 1.0; norZ[i] = 0.0;
        binX[i] = 0.0; binY[i] = 0.0; binZ[i] = 1.0;
        vNorX[i] = vNorY[i] = vNorZ[i] = 0.0;
        vBinX[i] = vBinY[i] = vBinZ[i] = 0.0;
      }
    #else // 2d stefan swimmer
      const std::array<Real ,6> curvature_points = {
          0, .15*length, .4*length, .65*length, .9*length, length
      };
      const std::array<Real ,6> curvature_values = {
          0.82014/length, 1.46515/length, 2.57136/length,
          3.75425/length, 5.09147/length, 5.70449/length
      };
      curvScheduler.transition(time,0,1,curvature_values,curvature_values);
      // query the schedulers for current values
      curvScheduler.gimmeValues(time, curvature_points, Nm, rS, rC, vC);
      // construct the curvature
      for(int i=0; i<Nm; i++) {
        const Real darg = 2.*M_PI;
        const Real arg  = 2.*M_PI*(time -rS[i]/length) + M_PI*phaseShift;
        rK[i] =   rC[i]*std::sin(arg);
        vK[i] =   vC[i]*std::sin(arg) + rC[i]*std::cos(arg)*darg;
      }

      // solve frenet to compute midline parameters
      IF2D_Frenet2D::solve(Nm, rS, rK, vK, rX, rY, vX, vY, norX, norY, vNorX, vNorY);
    #endif
  }
};

Naca::Naca(SimulationData&s, ArgumentParser&p) : Fish(s, p)
{
  Apitch = p("-Apitch").asDouble(0.0)*M_PI/180; //aplitude of sinusoidal pitch angle
  Fpitch = p("-Fpitch").asDouble(0.0)         ; //frequency
  Mpitch = p("-Mpitch").asDouble(0.0)*M_PI/180; //mean angle
  Fheave = p("-Fheave").asDouble(0.0)         ; //frequency of rowing motion
  Aheave = p("-Aheave").asDouble(0.0)*length  ; //amplitude (NON DIMENSIONAL)
  tAccel = p("-tAccel").asDouble(-1);
  fixedCenterDist = p("-fixedCenterDist").asDouble(0);
  const Real thickness = p("-tRatio").asDouble(0.12);
  myFish = new NacaMidlineData(length, sim.hmin, sim.extents[2], thickness);
  if( sim.rank == 0 && sim.verbose) printf("[CUP3D] - NacaData Nm=%d L=%f t=%f A=%f w=%f xvel=%f yvel=%f tAccel=%f fixedCenterDist=%f\n",myFish->Nm, (double)length, (double)thickness, (double)Apitch, (double)Fpitch, (double)transVel_imposed[0], (double)transVel_imposed[1], (double)tAccel, (double)fixedCenterDist);
  //only allow rotation around z-axis and translation in xy-plane
  bBlockRotation[0] = true;
  bBlockRotation[1] = true;
  bForcedInSimFrame[2] = true;
}

void Naca::computeVelocities()
{
  const Real omegaAngle = 2*M_PI*Fpitch;
  const Real angle = Mpitch + Apitch*std::sin(omegaAngle*sim.time);
  const Real omega = Apitch*omegaAngle*std::cos(omegaAngle*sim.time);
  // angular velocity
  angVel[0] = 0;
  angVel[1] = 0;
  angVel[2] = omega;

  // heaving motion
  const Real v_heave = -2.0*M_PI*Fheave*Aheave*std::sin(2*M_PI*Fheave*sim.time);
  if( sim.time < tAccel )
  {
    // linear velocity (due to rotation-axis != CoM)
    transVel[0] = (1.0-sim.time/tAccel)*0.01*transVel_imposed[0] + (sim.time/tAccel)*transVel_imposed[0] - fixedCenterDist*length*omega*std::sin(angle);
    transVel[1] = (1.0-sim.time/tAccel)*0.01*transVel_imposed[1] + (sim.time/tAccel)*transVel_imposed[1] + fixedCenterDist*length*omega*std::cos(angle) + v_heave;
    transVel[2] = 0.0;
  }
  else
  {
    // linear velocity (due to rotation-axis != CoM)
    transVel[0] = transVel_imposed[0] - fixedCenterDist*length*omega*std::sin(angle);
    transVel[1] = transVel_imposed[1] + fixedCenterDist*length*omega*std::cos(angle) + v_heave;
    transVel[2] = 0.0;
  }
}

using intersect_t = std::vector<std::vector<VolumeSegment_OBB*>>;
void Naca::writeSDFOnBlocks(std::vector<VolumeSegment_OBB> & vSegments)
{
  #pragma omp parallel
  {
    PutNacaOnBlocks putfish(myFish, position, quaternion);
    #pragma omp for
    for (size_t j=0 ; j < MyBlockIDs.size(); j++)
    {
      std::vector<VolumeSegment_OBB*> S;
      for (size_t k = 0 ; k < MySegments[j].size() ; k++)
      {
        VolumeSegment_OBB*const ptr  = & vSegments[MySegments[j][k]];
        S.push_back(ptr);
      }
      if(S.size() > 0)
      {
        ObstacleBlock*const block = obstacleBlocks[MyBlockIDs[j].blockID];
        putfish(MyBlockIDs[j].h,
                MyBlockIDs[j].origin_x,
                MyBlockIDs[j].origin_y,
                MyBlockIDs[j].origin_z, block, S);
      }
    }
  }
}

void Naca::updateLabVelocity( int nSum[3], Real uSum[3] )
{
  // heaving motion
  const Real v_heave = -2.0*M_PI*Fheave*Aheave*std::sin(2*M_PI*Fheave*sim.time);

  if(bFixFrameOfRef[0])
  {
   (nSum[0])++; 
   if( sim.time < tAccel ) uSum[0] -= (1.0-sim.time/tAccel)*0.01*transVel_imposed[0] + (sim.time/tAccel)*transVel_imposed[0];
   else                    uSum[0] -=                   transVel_imposed[0];
  }
  if(bFixFrameOfRef[1])
  {
   (nSum[1])++; 
   if( sim.time < tAccel )uSum[1] -= (1.0-sim.time/tAccel)*0.01*transVel_imposed[1] + (sim.time/tAccel)*transVel_imposed[1] + v_heave;
   else                   uSum[1] -=                   transVel_imposed[1] + v_heave; 
  }
  if(bFixFrameOfRef[2])
  {
   (nSum[2])++; 
   uSum[2] -= transVel[2]; 
  }

}

void Naca::update()
{
  const Real angle_2D =  Mpitch + Apitch * std::cos( 2*M_PI * (Fpitch*sim.time) );
  quaternion[0] = std::cos(0.5*angle_2D);
  quaternion[1] = 0;
  quaternion[2] = 0;
  quaternion[3] = std::sin(0.5*angle_2D);

  absPos[0] += sim.dt * transVel[0];
  absPos[1] += sim.dt * transVel[1];
  absPos[2] += sim.dt * transVel[2];

  position[0] += sim.dt * ( transVel[0] + sim.uinf[0] );
  // if user wants to keep airfoil in the mid plane then we just integrate
  // relative velocity (should be 0), otherwise we know that y velocity
  // is sinusoidal, therefore we can just use analytical form
  if(bFixFrameOfRef[1])
    position[1] += sim.dt * ( transVel[1] + sim.uinf[1] );
  else
    position[1] = sim.extents[1]/2 + Aheave * std::cos(2*M_PI*Fheave*sim.time);
  position[2] += sim.dt * ( transVel[2] + sim.uinf[2] );
}

CubismUP_3D_NAMESPACE_END
