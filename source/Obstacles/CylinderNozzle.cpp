//
//  Cubism3D
//  Copyright (c) 2023 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//

#include "CylinderNozzle.h"
#include "extra/ObstacleLibrary.h"

#include <Cubism/ArgumentParser.h>

CubismUP_3D_NAMESPACE_BEGIN
using namespace cubism;

CylinderNozzle::CylinderNozzle(SimulationData&s, ArgumentParser &p): 
  Cylinder(s, p),
  Nactuators     ( p("-Nactuators"    ).asInt(2)               ),
  actuator_theta ( p("-actuator_theta").asDouble(10.)*M_PI/180.),
  regularizer    ( p("-regularizer"   ).asDouble(1.0)          ),
  ccoef( p("-ccoef").asDouble(0.1) )
{
    actuators.resize(Nactuators,0.);
    actuatorSchedulers.resize(Nactuators);
    actuators_prev_value.resize(Nactuators);
    actuators_next_value.resize(Nactuators);
}


void CylinderNozzle::finalize()
{
  const double cd = force[0] / (0.5*transVel[0]*transVel[0]*2*radius*2*halflength);
  fx_integral += -std::fabs(cd)*sim.dt;

  const Real transition_duration = 0.1;
  for (size_t idx = 0 ; idx < actuators.size(); idx++)
  {
    Real dummy;
    actuatorSchedulers[idx].transition (sim.time,t_change,t_change+transition_duration,actuators_prev_value[idx],actuators_next_value[idx]);
    actuatorSchedulers[idx].gimmeValues(sim.time,actuators[idx],dummy);
  }

  const auto & vInfo = sim.vel->getBlocksInfo();
  const Real dtheta = 2*M_PI/Nactuators;
  const Real Cx = position[0];
  const Real Cy = position[1];
  const Real Cz = position[2];
  const Real Uact_max = ccoef * pow(transVel[0]*transVel[0] + transVel[1]*transVel[1] + transVel[2]*transVel[2],0.5);
  for(size_t i=0; i<vInfo.size(); i++)
  {
    const auto & info = vInfo[i];
    if(obstacleBlocks[info.blockID] == nullptr) continue; //obst not in block
    ObstacleBlock& o = * obstacleBlocks[info.blockID];
    auto & __restrict__ UDEF = o.udef;

    for(int iz=0; iz<ScalarBlock::sizeZ; iz++)
    for(int iy=0; iy<ScalarBlock::sizeY; iy++)
    for(int ix=0; ix<ScalarBlock::sizeX; ix++)
    {
      UDEF[iz][iy][ix][0] = 0.0;
      UDEF[iz][iy][ix][1] = 0.0;
      UDEF[iz][iy][ix][2] = 0.0;
      Real p[3];
      info.pos(p, ix, iy, iz);
      const Real x = p[0]-Cx;
      const Real y = p[1]-Cy;
      const Real z = p[2]-Cz;
      const Real r = x*x+y*y;
      if (r > (radius+2*info.h)*(radius+2*info.h) || r < (radius-2*info.h)*(radius-2*info.h)) continue;
      if (std::fabs(z) > 0.99*halflength) continue;
      //if (std::fabs(z) > 0.75*halflength) continue;

      Real theta = atan2(y,x);
      if (theta < 0) theta += 2.*M_PI;

      int idx = round(theta / dtheta); //this is the closest actuator
      if (idx == Nactuators) idx = 0;  //periodic around the cylinder

      const Real theta0 = idx * dtheta;

      const Real phi = theta - theta0;
      if ( std::fabs(phi) < 0.5*actuator_theta || (idx == 0 && std::fabs(phi-2*M_PI) < 0.5*actuator_theta))
      {
        const Real rr = radius / pow(r,0.5);
        const Real ur = Uact_max*rr*actuators[idx]*cos(M_PI*phi/actuator_theta);
        UDEF[iz][iy][ix][0] = ur * cos(theta);
        UDEF[iz][iy][ix][1] = ur * sin(theta);
        UDEF[iz][iy][ix][2] = 0.0;
      }
    }
  }
}

void CylinderNozzle::act( std::vector<Real> action, const int agentID)
{
   t_change = sim.time;
   if(action.size() != actuators.size())
   {
     std::cerr << "action size needs to be equal to actuators\n";
     fflush(0);
     abort();
   }

   bool bounded = false;
   while (bounded == false)
   {
       bounded = true;
       Real Q = 0;
       for (size_t i = 0 ; i < action.size() ; i ++)
       {
            Q += action[i];
       }
       Q /= action.size();
       for (size_t i = 0 ; i < action.size() ; i ++)
       {
            action[i] -= Q;
            if (std::fabs(action[i]) > 1.0) bounded = false;
            action[i] = std::max(action[i],-1.0);
            action[i] = std::min(action[i],+1.0);
       }
   }

   for (size_t i = 0 ; i < action.size() ; i ++)
   {
     actuators_prev_value[i] = actuators[i];
     actuators_next_value[i] = action   [i];
   }
}


Real CylinderNozzle::reward(const int agentID)
{
  Real retval = fx_integral / 0.1; //0.1 is the action times
  fx_integral = 0;
  Real regularizer_sum = 0.0;
  for (size_t idx = 0 ; idx < actuators.size(); idx++)
  {
    regularizer_sum += actuators[idx]*actuators[idx];
  }
  regularizer_sum = pow(regularizer_sum,0.5)/actuators.size(); //O(1)
  const double c = -regularizer;
  return retval + c*regularizer_sum;
}

std::vector<Real> CylinderNozzle::state(const int agentID)
{
  std::vector<Real> S;
  const int bins = 16;
  const Real bins_theta = 10*M_PI/180.0;
  const Real dtheta = 2.*M_PI / bins;
  std::vector<int>  n_s(bins,0.0);
  std::vector<Real> p_s(bins,0.0);
  std::vector<Real> o_s(bins,0.0);
  for(auto & block : obstacleBlocks) if(block not_eq nullptr)
  {
    for(int i=0; i<block->nPoints; i++)
    {
      const Real x = block->pX[i] - position[0];
      const Real y = block->pY[i] - position[1];
      const Real z = block->pZ[i] - position[2];
      if (std::fabs(z) > 0.75*halflength) continue;

      Real theta = atan2(y,x);
      if (theta < 0) theta += 2.*M_PI;
      int idx = round(theta / dtheta); //this is the closest actuator
      if (idx == bins) idx = 0;  //periodic around the cylinder
      const Real theta0 = idx * dtheta;
      const Real phi = theta - theta0;
      if ( std::fabs(phi) < 0.5*bins_theta || (idx == 0 && std::fabs(phi-2*M_PI) < 0.5*bins_theta))
      {
        const Real p = block->P[i];
        const Real om = block->omegaZ[i];
        n_s[idx] ++;
        p_s[idx] += p;
        o_s[idx] += om;
      }
    }
  }

  MPI_Allreduce(MPI_IN_PLACE,n_s.data(),n_s.size(),MPI_INT ,MPI_SUM,sim.comm);
  for (int idx = 0 ; idx < bins; idx++)
  {
    if (n_s[idx] == 0) continue;
    p_s[idx] /= n_s[idx];
    o_s[idx] /= n_s[idx];
  }

  for (int idx = 0 ; idx < bins; idx++) S.push_back(p_s[idx]);
  for (int idx = 0 ; idx < bins; idx++) S.push_back(o_s[idx]);
  MPI_Allreduce(MPI_IN_PLACE,  S.data(),  S.size(),MPI_Real,MPI_SUM,sim.comm);
  S.push_back(-force[0]/(2*halflength));
  S.push_back(-force[1]/(2*halflength));
  const Real Re = std::fabs(transVel[0])*(2*radius)/sim.nu;
  S.push_back(Re);
  S.push_back(ccoef);

  if (sim.rank ==0 )
    for (size_t i = 0 ; i < S.size() ; i++) std::cout << S[i] << " ";

  return S;
}

CubismUP_3D_NAMESPACE_END
