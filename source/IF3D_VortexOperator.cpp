//
//  CubismUP_3D
//
//  Written by Guido Novati ( novatig@ethz.ch ).
//  This file started as an extension of code written by Wim van Rees
//  Copyright (c) 2017 ETHZ. All rights reserved.
//

#include "IF3D_VortexOperator.h"
#include <chrono>
IF3D_VortexOperator::IF3D_VortexOperator(FluidGridMPI * g, ArgumentParser & p, const Real*const u) : IF3D_ObstacleOperator(g, p, u), created(false)
{
  v_max = p("-vmax").asDouble();
}

void IF3D_VortexOperator::create(const int step_id,const double time, const double dt, const Real *Uinf)
{
  if(time<8. || created) return;
  printf("Crating a vortex in %f %f %f core size %f max vel %f\n", position[0],position[1],position[2],length,v_max);
  #pragma omp parallel
  {
    const int N = vInfo.size();
    const Real alpha=1.25643, fac=1.39795, eps=numeric_limits<Real>::epsilon();
    #pragma omp for schedule(static)
    for(int i=0; i<vInfo.size(); i++) {
      BlockInfo info = vInfo[i];
      FluidBlock& b = *(FluidBlock*)info.ptrBlock;
      for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
      for(int iy=0; iy<FluidBlock::sizeY; ++iy)
      for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
        double p[3];
        info.pos(p, ix, iy, iz);
        p[0]-=position[0];
        p[1]-=position[1];
        p[2]-=position[2];
        const double theta = std::atan2(p[1], p[0]);
        const double r = std::sqrt(p[0]*p[0] + p[1]*p[1])/length;
        //if(r>3) continue;
        const double z = std::fabs(p[2])/0.05;
        //if(z>3) continue;
        const double arg = std::max((double)-20, -alpha*r*r -z*z);
        const double vTheta = 0.1*fac/(r+eps)*(1-std::exp(arg));
        b(ix,iy,iz).u -= vTheta*std::sin(theta);
        b(ix,iy,iz).v += vTheta*std::cos(theta);
      }
    }
  }
  created = true;
}

void IF3D_VortexOperator::update(const int stepID, const double t, const double dt, const Real *Uinf)
{}

void IF3D_VortexOperator::computeVelocities(const Real* Uinf) {}
void IF3D_VortexOperator::computeForces(const int stepID, const double time, const double dt, const Real* Uinf, const double NU, const bool bDump) {}
void IF3D_VortexOperator::Accept(ObstacleVisitor * visitor) {}
void IF3D_VortexOperator::finalize(const int step_id,const double time, const double dt, const Real *Uinf) {}
void IF3D_VortexOperator::execute(const int iAgent, const double time, const vector<double>a) {}
void IF3D_VortexOperator::interpolateOnSkin(const double time, const int stepID, bool dumpWake) {}
void IF3D_VortexOperator::getSkinsAndPOV(Real& x, Real& y, Real& th, Real*& pXL, Real*& pYL, Real*& pXU, Real*& pYU, int& Npts)
{
  pXL = pYL = pXU = pYU = nullptr;
  Npts = 0;
}
void IF3D_VortexOperator::computeDiagnostics(const int stepID, const double time, const Real* Uinf, const double lambda) {}
