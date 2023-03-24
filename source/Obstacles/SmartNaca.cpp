//
//  CubismUP_3D
//  Copyright (c) 2023 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//

#include "SmartNaca.h"
#include "FishLibrary.h"

using namespace cubism;

CubismUP_3D_NAMESPACE_BEGIN

SmartNaca::SmartNaca(SimulationData&s, ArgumentParser&p): Naca(s,p), Nactuators ( p("-Nactuators").asInt(2)),actuator_ds( p("-actuatords").asDouble(0.05)),thickness(p("-tRatio").asDouble(0.12))
{
  actuators.resize(Nactuators,0.);
  actuatorSchedulers.resize(Nactuators);
  actuators_prev_value.resize(Nactuators);
  actuators_next_value.resize(Nactuators);
}

void SmartNaca::finalize()
{
  const Real * const rS = myFish->rS;
  /*const*/ Real * const rX = myFish->rX;
  /*const*/ Real * const rY = myFish->rY;
  /*const*/ Real * const rZ = myFish->rZ;
  /*const*/ Real * const norX = myFish->norX;
  /*const*/ Real * const norY = myFish->norY;
  /*const*/ Real * const norZ = myFish->norZ;
  /*const*/ Real * const binX = myFish->binX;
  /*const*/ Real * const binY = myFish->binY;
  /*const*/ Real * const binZ = myFish->binZ;
  const Real * const width = myFish->width;
  //const Real * const height = myFish->height;

  PutFishOnBlocks putfish(myFish, position, quaternion);
  for (int ss = 0 ; ss < myFish->Nm; ss++)
  {
    Real x[3] = {rX[ss],rY[ss],rZ[ss]};
    Real n[3] = {norX[ss],norY[ss],norZ[ss]};
    Real b[3] = {binX[ss],binY[ss],binZ[ss]};
    putfish.changeToComputationalFrame(x);
    putfish.changeVelocityToComputationalFrame(n);
    putfish.changeVelocityToComputationalFrame(b);
    rX[ss] = x[0]; 
    rY[ss] = x[1]; 
    rZ[ss] = x[2]; 
    norX[ss] = n[0]; 
    norY[ss] = n[1]; 
    norZ[ss] = n[2]; 
    binX[ss] = b[0]; 
    binY[ss] = b[1]; 
    binZ[ss] = b[2]; 
  }

  #if 1
  //dummy actuator values for testing
  static bool visited = false;
  #if 1
  if (sim.time > 2.0 && visited == false)
  {
    visited = true;
    std::vector<Real> q(actuators.size());
    for (int i = 0 ; i < (int)actuators.size(); i ++) q[i] = 0.25*(2*(i+1)%2-1);
    act(q,0);
  }
  #endif

  const Real transition_duration = 0.1;
  Real tot = 0.0;
  for (size_t idx = 0 ; idx < actuators.size(); idx++)
  {
    Real dummy;
    actuatorSchedulers[idx].transition (sim.time,t_change,t_change+transition_duration,actuators_prev_value[idx],actuators_next_value[idx]);
    actuatorSchedulers[idx].gimmeValues(sim.time,actuators[idx],dummy);
    tot += std::fabs(actuators[idx]);
  }
  const double cd = force[0] / (0.5*transVel[0]*transVel[0]*thickness*thickness);
  fx_integral += -std::fabs(cd)*sim.dt;
  if (tot < 1e-21) return;


  std::vector<int>       ix_store;
  std::vector<int>       iy_store;
  std::vector<int>       iz_store;
  std::vector<long long> id_store;
  std::vector<Real>      nx_store;
  std::vector<Real>      ny_store;
  std::vector<Real>      nz_store;

  for (const auto & info : sim.vel->getBlocksInfo())
  {
    if(obstacleBlocks[info.blockID] == nullptr) continue; //obst not in block
    ObstacleBlock& o = * obstacleBlocks[info.blockID];
    auto & __restrict__ UDEF = o.udef;
    auto & __restrict__ SDF  = o.sdfLab;

    for(int iz=0; iz<ScalarBlock::sizeZ; iz++)
    for(int iy=0; iy<ScalarBlock::sizeY; iy++)
    for(int ix=0; ix<ScalarBlock::sizeX; ix++)
    {
      if ( SDF[iz+1][iy+1][ix+1] > info.h || SDF[iz+1][iy+1][ix+1] < -info.h) continue;
      UDEF[iz][iy][ix][0] = 0.0;
      UDEF[iz][iy][ix][1] = 0.0;
      UDEF[iz][iy][ix][2] = 0.0;
      Real p[3];
      info.pos(p, ix, iy, iz);

      if ( std::fabs(p[2]-position[2])> 0.40*thickness ) continue;

      //find closest surface point to analytical expression
      int  ss_min = 0;
      int  sign_min = 0;
      Real dist_min = 1e10;

      for (int ss = 0 ; ss < myFish->Nm; ss++)
      {
          Real Pp [3] = {rX[ss]+width[ss]*norX[ss],rY[ss]+width[ss]*norY[ss],0};
          Real Pm [3] = {rX[ss]-width[ss]*norX[ss],rY[ss]-width[ss]*norY[ss],0};
          const Real dp = pow(Pp[0]-p[0],2)+pow(Pp[1]-p[1],2);
          const Real dm = pow(Pm[0]-p[0],2)+pow(Pm[1]-p[1],2);
          if (dp < dist_min)
          {
            sign_min = 1;
            dist_min = dp;
            ss_min = ss;
          }
          if (dm < dist_min)
          {
            sign_min = -1;
            dist_min = dm;
            ss_min = ss;
          }
      }

      const Real smax = rS[myFish->Nm-1]-rS[0];
      const Real ds   = 2*smax/Nactuators;
      const Real current_s = rS[ss_min];
      if (current_s < 0.01*length || current_s > 0.99*length) continue;
      int idx = (current_s / ds); //this is the closest actuator
      const Real s0 = 0.5*ds + idx * ds;
      if (sign_min == -1) idx += Nactuators/2;

      if (std::fabs( current_s - s0 ) < 0.5*actuator_ds*length)
      {
        Real nx = 0.5*(SDF[iz+1][iy+1][ix+2]-SDF[iz+1][iy+1][ix+0]);
        Real ny = 0.5*(SDF[iz+1][iy+2][ix+1]-SDF[iz+1][iy+0][ix+1]);
        Real nz = 0.5*(SDF[iz+2][iy+1][ix+1]-SDF[iz+0][iy+1][ix+1]);
        const Real nn = pow(nx*nx+ny*ny+nz*nz+1e-21,-0.5);
        nx *= nn;
        ny *= nn;
        nz *= nn;
        const double c = 1.0 - std::fabs(current_s - s0)/ (0.5*actuator_ds*length);
        UDEF[iz][iy][ix][0] = c*actuators[idx]*nx;
        UDEF[iz][iy][ix][1] = c*actuators[idx]*ny;
        UDEF[iz][iy][ix][2] = c*actuators[idx]*nz;
        ix_store.push_back(ix);
        iy_store.push_back(iy);
        iz_store.push_back(iz);
        id_store.push_back(info.blockID);
        nx_store.push_back(nx);
        ny_store.push_back(ny);
        nz_store.push_back(nz);
      }
    }
  }

  //Compute surface and total mass flux
  Real surface   = 0.0;
  Real mass_flux = 0.0;
  #pragma omp parallel for reduction(+: surface,mass_flux)
  for (const auto & info : sim.vel->getBlocksInfo())
  {
    if(obstacleBlocks[info.blockID] == nullptr) continue; //obst not in block
    ObstacleBlock& o = * obstacleBlocks[info.blockID];
    auto & __restrict__ UDEF = o.udef;
    auto & __restrict__ SDF  = o.sdfLab;
    auto & __restrict__ CHI  = o.chi;
    const Real invh = 1.0/info.h;
    const Real h3 = info.h*info.h*info.h;

    for(int iz=0; iz<ScalarBlock::sizeZ; iz++)
    for(int iy=0; iy<ScalarBlock::sizeY; iy++)
    for(int ix=0; ix<ScalarBlock::sizeX; ix++)
    {
      Real nx = 0.5*(SDF[iz+1][iy+1][ix+2]-SDF[iz+1][iy+1][ix+0]);
      Real ny = 0.5*(SDF[iz+1][iy+2][ix+1]-SDF[iz+1][iy+0][ix+1]);
      Real nz = 0.5*(SDF[iz+2][iy+1][ix+1]-SDF[iz+0][iy+1][ix+1]);
      const Real nn = pow(nx*nx+ny*ny+nz*nz+1e-21,-0.5);
      nx *= nn;
      ny *= nn;
      nz *= nn;
      const Real dchidx = invh * ( (ix > 0 && ix < ScalarBlock::sizeX-1) ? 0.5*(CHI[iz][iy][ix+1]-CHI[iz][iy][ix-1]):(ix == 0 ? (-1.5*CHI[iz][iy][ix]+2.0*CHI[iz][iy][ix+1]-0.5*CHI[iz][iy][ix+2]) : (1.5*CHI[iz][iy][ix]-2.0*CHI[iz][iy][ix-1]+0.5*CHI[iz][iy][ix-2]) ) );
      const Real dchidy = invh * ( (iy > 0 && iy < ScalarBlock::sizeY-1) ? 0.5*(CHI[iz][iy+1][ix]-CHI[iz][iy-1][ix]):(iy == 0 ? (-1.5*CHI[iz][iy][ix]+2.0*CHI[iz][iy+1][ix]-0.5*CHI[iz][iy+2][ix]) : (1.5*CHI[iz][iy][ix]-2.0*CHI[iz][iy-1][ix]+0.5*CHI[iz][iy-2][ix]) ) );
      const Real dchidz = invh * ( (iz > 0 && iz < ScalarBlock::sizeZ-1) ? 0.5*(CHI[iz+1][iy][ix]-CHI[iz-1][iy][ix]):(iz == 0 ? (-1.5*CHI[iz][iy][ix]+2.0*CHI[iz+1][iy][ix]-0.5*CHI[iz+2][iy][ix]) : (1.5*CHI[iz][iy][ix]-2.0*CHI[iz-1][iy][ix]+0.5*CHI[iz-2][iy][ix]) ) );
      const Real fac = (dchidx*nx+dchidy*ny+dchidz*nz)*h3;
      mass_flux += fac*(UDEF[iz][iy][ix][0]*nx+UDEF[iz][iy][ix][1]*ny+UDEF[iz][iy][ix][2]*nz);
      surface   += fac;
    }
  }
  Real Qtot [2] = {mass_flux,surface};
  MPI_Allreduce(MPI_IN_PLACE,Qtot,2,MPI_Real,MPI_SUM,sim.comm);
  const Real uMean = Qtot[0]/Qtot[1];
  //if (sim.rank == 0) std::cout << "Q=" << Qtot[0] << " " <<  Qtot[1] << std::endl;

  //Substract total mass flux (divided by surface) from actuator velocities
  #pragma omp parallel for
  for (size_t idx = 0 ; idx < id_store.size(); idx++)
  {
    const long long blockID = id_store[idx];
    const int ix            = ix_store[idx];
    const int iy            = iy_store[idx];
    const int iz            = iz_store[idx];
    const Real nx           = nx_store[idx];
    const Real ny           = ny_store[idx];
    const Real nz           = nz_store[idx];
    ObstacleBlock& o = * obstacleBlocks[blockID];
    auto & __restrict__ UDEF = o.udef;
    UDEF[iz][iy][ix][0] -= uMean*nx;
    UDEF[iz][iy][ix][1] -= uMean*ny;
    UDEF[iz][iy][ix][2] -= uMean*nz;
  }
  #endif
}

void SmartNaca::act( std::vector<Real> action, const int agentID)
{
  t_change = sim.time;
  if(action.size() != actuators.size())
  {
    std::cerr << "action size needs to be equal to actuators\n";
    fflush(0);
    abort();
  }
  for (size_t i = 0 ; i < action.size() ; i ++)
  {
    actuators_prev_value[i] = actuators[i];
    actuators_next_value[i] = action   [i];
  }
}

Real SmartNaca::reward(const int agentID)
{
  Real retval = fx_integral / 0.1; //0.1 is the action times
  fx_integral = 0;
  Real regularizer = 0.0;
  for (size_t idx = 0 ; idx < actuators.size(); idx++)
  {
    regularizer += actuators[idx]*actuators[idx];
  }
  regularizer = pow(regularizer,0.5)/actuators.size();
  return retval - 0.1*regularizer;
}

std::vector<Real> SmartNaca::state(const int agentID)
{
  std::vector<Real> S;
  #if 0
  const int bins = 64;

  const Real dtheta = 2.*M_PI / bins;
  std::vector<int>   n_s   (bins,0.0);
  std::vector<Real>  p_s   (bins,0.0);
  std::vector<Real> fX_s   (bins,0.0);
  std::vector<Real> fY_s   (bins,0.0);
  for(auto & block : obstacleBlocks) if(block not_eq nullptr)
  {
    for(size_t i=0; i<block->n_surfPoints; i++)
    {
      const Real x     = block->x_s[i] - origC[0];
      const Real y     = block->y_s[i] - origC[1];
      const Real ang   = atan2(y,x);
      const Real theta = ang < 0 ? ang + 2*M_PI : ang;
      const Real p     = block->p_s[i];
      const Real fx    = block->fX_s[i];
      const Real fy    = block->fY_s[i];
      const int idx = theta / dtheta;
      n_s [idx] ++;
      p_s [idx] += p;
      fX_s[idx] += fx;
      fY_s[idx] += fy;
    }
  }

  MPI_Allreduce(MPI_IN_PLACE,n_s.data(),n_s.size(),MPI_INT ,MPI_SUM,sim.comm);
  for (int idx = 0 ; idx < bins; idx++)
  {
    p_s [idx] /= n_s[idx];
    fX_s[idx] /= n_s[idx];
    fY_s[idx] /= n_s[idx];
  }

  for (int idx = 0 ; idx < bins; idx++) S.push_back( p_s[idx]);
  for (int idx = 0 ; idx < bins; idx++) S.push_back(fX_s[idx]);
  for (int idx = 0 ; idx < bins; idx++) S.push_back(fY_s[idx]);
  MPI_Allreduce(MPI_IN_PLACE,  S.data(),  S.size(),MPI_Real,MPI_SUM,sim.comm);
  S.push_back(forcex);
  S.push_back(forcey);
  S.push_back(torque);

  if (sim.rank ==0 )
    for (size_t i = 0 ; i < S.size() ; i++) std::cout << S[i] << " ";

  #endif
  return S;
}

CubismUP_3D_NAMESPACE_END
