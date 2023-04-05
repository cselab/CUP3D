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
  regularizer    ( p("-regularizer"   ).asDouble(1.0)          )
{
    actuators.resize(Nactuators,0.);
    actuatorSchedulers.resize(Nactuators);
    actuators_prev_value.resize(Nactuators);
    actuators_next_value.resize(Nactuators);
    #if 0
    if (sim.rank == 0)
    {
      std::string line;
      ifstream myfile;
      myfile.open ("actions0.txt",std::ios_base::in);

      while (std::getline(myfile, line))
      {
          std::istringstream iss(line);
          std::cout << "--> ";
          double temp;
          iss >> temp;
          t_action_taken.push_back(temp-100);
          std::cout << temp << " ";
          for( size_t a = 0; a<actuators.size(); a++ )
          {
            iss >> temp;
            std::cout << temp << " ";
            action_taken.push_back(temp);
          }
          std::cout << std::endl;
      }
      myfile.close();
    }
    int n = (int) t_action_taken.size();
    MPI_Bcast(&n, 1, MPI_INT, 0, sim.comm );
    t_action_taken.resize(n);
    action_taken.resize(actuators.size()*n);
    MPI_Bcast(t_action_taken.data(), t_action_taken.size(), MPI_DOUBLE, 0, sim.comm );
    MPI_Bcast(  action_taken.data(),   action_taken.size(), MPI_DOUBLE, 0, sim.comm );
    #endif
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

//  if (sim.time >= t_action_taken[curr_idx])
//  {
//    std::vector<Real> a (actuators.size());
//    for (size_t i = 0 ; i < actuators.size(); i++)
//      a[i] = action_taken[curr_idx*actuators.size()+i];
//    act(a,0);
//    curr_idx++;
//  }

  const auto & vInfo = sim.vel->getBlocksInfo();
  const Real dtheta = 2*M_PI/Nactuators;
  const Real Cx = position[0];
  const Real Cy = position[1];
  const Real Cz = position[2];
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
      if (std::fabs(z) > 0.4*halflength) continue;

      Real theta = atan2(y,x);
      if (theta < 0) theta += 2.*M_PI;

      int idx = round(theta / dtheta); //this is the closest actuator
      if (idx == Nactuators) idx = 0;  //periodic around the cylinder

      const Real theta0 = idx * dtheta;

      const Real phi = theta - theta0;
      if ( std::fabs(phi) < 0.5*actuator_theta || (idx == 0 && std::fabs(phi-2*M_PI) < 0.5*actuator_theta))
      {
        const Real rr = pow(r,0.5);
        //const Real ur = 0.01*actuators[idx]/rr*cos(M_PI*phi/actuator_theta);
        const Real ur = 0.005*actuators[idx]/rr*cos(M_PI*phi/actuator_theta);
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
  Real Q = 0;
  for (size_t i = 0 ; i < action.size() ; i ++)
  {
    actuators_prev_value[i] = actuators[i];
    actuators_next_value[i] = action   [i];
    Q += action[i];
  }
  Q /= actuators.size();
  for (size_t i = 0 ; i < action.size() ; i ++) actuators_next_value[i] -= Q;
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
  std::vector<int>   n_s   (bins,0.0);
  std::vector<Real>  p_s   (bins,0.0);
  std::vector<Real> fX_s   (bins,0.0);
  std::vector<Real> fY_s   (bins,0.0);
  for(auto & block : obstacleBlocks) if(block not_eq nullptr)
  {
    for(int i=0; i<block->nPoints; i++)
    {
      const Real x = block->pX[i] - position[0];
      const Real y = block->pY[i] - position[1];
      Real theta = atan2(y,x);
      if (theta < 0) theta += 2.*M_PI;
      int idx = round(theta / dtheta); //this is the closest actuator
      if (idx == bins) idx = 0;  //periodic around the cylinder
      const Real theta0 = idx * dtheta;
      const Real phi = theta - theta0;
      if ( std::fabs(phi) < 0.5*bins_theta || (idx == 0 && std::fabs(phi-2*M_PI) < 0.5*bins_theta))
      {
        const Real p     = block->P[i];
        const Real fx    = block->fX[i];
        const Real fy    = block->fY[i];
        n_s [idx] ++;
        p_s [idx] += p;
        fX_s[idx] += fx;
        fY_s[idx] += fy;
      }
    }
  }

  MPI_Allreduce(MPI_IN_PLACE,n_s.data(),n_s.size(),MPI_INT ,MPI_SUM,sim.comm);
  for (int idx = 0 ; idx < bins; idx++)
  {
    if (n_s[idx] == 0) continue;
    p_s [idx] /= n_s[idx];
    fX_s[idx] /= n_s[idx];
    fY_s[idx] /= n_s[idx];
  }

  for (int idx = 0 ; idx < bins; idx++) S.push_back( p_s[idx]);
  for (int idx = 0 ; idx < bins; idx++) S.push_back(fX_s[idx]);
  for (int idx = 0 ; idx < bins; idx++) S.push_back(fY_s[idx]);
  MPI_Allreduce(MPI_IN_PLACE,  S.data(),  S.size(),MPI_Real,MPI_SUM,sim.comm);
  S.push_back(force[0]);
  S.push_back(force[1]);

  if (sim.rank ==0 )
    for (size_t i = 0 ; i < S.size() ; i++) std::cout << S[i] << " ";

  return S;
}



CubismUP_3D_NAMESPACE_END
