//
//  CubismUP_3D
//
//  Written by Guido Novati ( novatig@ethz.ch ).
//  Copyright (c) 2017 ETHZ. All rights reserved.
//

#include "IF3D_ObstacleOperator.h"
#include "GenericOperator.h"
#include "BufferedLogger.h"

void IF3D_ObstacleOperator::_computeUdefMoments(double lin_momenta[3],
                                        double ang_momenta[3], const double CoM[3])
{
  const double h   = vInfo[0].h_gridpoint;
  const double dv  = std::pow(vInfo[0].h_gridpoint, 3);
  constexpr double eps = std::numeric_limits<Real>::epsilon();
  (void)eps;

  { //sum linear momenta to figure out velocity and mass
    double V=0, lm0=0, lm1=0, lm2=0; //linear momenta
    #pragma omp parallel for schedule(dynamic) reduction(+:V,lm0,lm1,lm2)
    for(size_t i=0; i<vInfo.size(); i++) {
      BlockInfo info = vInfo[i];
      const auto pos = obstacleBlocks.find(info.blockID);
      if(pos == obstacleBlocks.end()) continue;

      for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
      for(int iy=0; iy<FluidBlock::sizeY; ++iy)
      for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
        const double Xs = pos->second->chi[iz][iy][ix];
        if (Xs <= 0) continue;
        V   += Xs;
        lm0 += Xs * pos->second->udef[iz][iy][ix][0];
        lm1 += Xs * pos->second->udef[iz][iy][ix][1];
        lm2 += Xs * pos->second->udef[iz][iy][ix][2];
      }
    }

    double globals[4] = {0,0,0,0};
    double locals[4] = {lm0,lm1,lm2,V};
    MPI_Allreduce(locals, globals, 4, MPI_DOUBLE, MPI_SUM, grid->getCartComm());

    assert(globals[3] > std::numeric_limits<Real>::epsilon());
    #ifdef _VERBOSE_
      if (!rank)
      printf("Discrepancy in computed volume during correction = %g.\n",
        std::fabs(computed_vol- globals[3]*dv));
    #endif
    lin_momenta[0] = globals[0]/globals[3];
    lin_momenta[1] = globals[1]/globals[3];
    lin_momenta[2] = globals[2]/globals[3];
    volume         = globals[3] * dv;
  }

  //sum angular momenta to figure out ang velocity and moments
  double J0=0, J1=0, J2=0, J3=0, J4=0, J5=0, am0=0, am1=0, am2=0;
  #pragma omp parallel for schedule(dynamic) reduction(+:J0,J1,J2,J3,J4,J5,am0,am1,am2)
  for(size_t i=0; i<vInfo.size(); i++) {
    BlockInfo info = vInfo[i];
    const auto pos = obstacleBlocks.find(info.blockID);
    if(pos == obstacleBlocks.end()) continue;

    for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
    for(int iy=0; iy<FluidBlock::sizeY; ++iy)
    for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
      const double Xs = pos->second->chi[iz][iy][ix];
      if (Xs <= 0) continue;
      assert(Xs>0);
      double p[3];
      info.pos(p, ix, iy, iz);
      p[0]-=CoM[0]; p[1]-=CoM[1]; p[2]-=CoM[2];
      const double u_ = pos->second->udef[iz][iy][ix][0];
      const double v_ = pos->second->udef[iz][iy][ix][1];
      const double w_ = pos->second->udef[iz][iy][ix][2];
      am0 += Xs * (p[1]*(w_-lin_momenta[2]) - p[2]*(v_-lin_momenta[1]));
      am1 += Xs * (p[2]*(u_-lin_momenta[0]) - p[0]*(w_-lin_momenta[2]));
      am2 += Xs * (p[0]*(v_-lin_momenta[1]) - p[1]*(u_-lin_momenta[0]));
      J0  += Xs * (p[1]*p[1]+p[2]*p[2]);
      J1  += Xs * (p[0]*p[0]+p[2]*p[2]);
      J2  += Xs * (p[0]*p[0]+p[1]*p[1]);
      J3  -= Xs *  p[0]*p[1];
      J4  -= Xs *  p[0]*p[2];
      J5  -= Xs *  p[1]*p[2];
    }
  }

  double globals[9] = {0,0,0,0,0,0,0,0,0};
  double locals[9] = {am0,am1,am2,J0,J1,J2,J3,J4,J5};
  MPI_Allreduce(locals, globals, 9, MPI_DOUBLE, MPI_SUM, grid->getCartComm());

  //if(bFixToPlanar)
  //{
  //  ang_momenta[0] = ang_momenta[1] = 0.0;
  //  ang_momenta[2] = globals[2]/globals[5]; // av2/j2
  //}
  //else
  {
    //solve avel = invJ \dot angMomentum, do not multiply by h^3, but by h for numerics
    const double AM[3] = {globals[0]*h, globals[1]*h, globals[2]*h};
    const double J_[6] = {globals[3]*h, globals[4]*h, globals[5]*h,
                        globals[6]*h, globals[7]*h, globals[8]*h};

    const double detJ = J_[0]*(J_[1]*J_[2] - J_[5]*J_[5])+
                      J_[3]*(J_[4]*J_[5] - J_[2]*J_[3])+
                      J_[4]*(J_[3]*J_[5] - J_[1]*J_[4]);
    const double invDetJ = 1./detJ;
    assert(std::fabs(detJ)>eps);
    const double invJ[6] = {
      invDetJ * (J_[1]*J_[2] - J_[5]*J_[5]),
      invDetJ * (J_[0]*J_[2] - J_[4]*J_[4]),
      invDetJ * (J_[0]*J_[1] - J_[3]*J_[3]),
      invDetJ * (J_[4]*J_[5] - J_[2]*J_[3]),
      invDetJ * (J_[3]*J_[5] - J_[1]*J_[4]),
      invDetJ * (J_[3]*J_[4] - J_[0]*J_[5])
    };

    ang_momenta[0] = invJ[0]*AM[0] + invJ[3]*AM[1] + invJ[4]*AM[2];
    ang_momenta[1] = invJ[3]*AM[0] + invJ[1]*AM[1] + invJ[5]*AM[2];
    ang_momenta[2] = invJ[4]*AM[0] + invJ[5]*AM[1] + invJ[2]*AM[2];
  }
  J[0] = globals[3] * dv;
  J[1] = globals[4] * dv;
  J[2] = globals[5] * dv;
  J[3] = globals[6] * dv;
  J[4] = globals[7] * dv;
  J[5] = globals[8] * dv;
  //assert(!errors);
  #ifdef _VERBOSE_
    const Real errors = std::max(std::fabs(globals[3]*dv-J[0]),
                          std::max(std::fabs(globals[4]*dv-J[1]),
                            std::max(std::fabs(globals[5]*dv-J[2]),
                              std::max(std::fabs(globals[6]*dv-J[3]),
                                std::max(std::fabs(globals[7]*dv-J[4]),
                                         std::fabs(globals[8]*dv-J[5]))))));
    if (!rank)
    printf("Max error in computed momenta during correction = %g.\n", errors);
  #endif
}

void IF3D_ObstacleOperator::_makeDefVelocitiesMomentumFree(const double CoM[3])
{
  _computeUdefMoments(transVel_correction, angVel_correction, CoM);
  #ifdef _VERBOSE_
   if(rank==0)
      printf("Correction of: lin mom [%f %f %f] ang mom [%f %f %f]\n",
        transVel_correction[0], transVel_correction[1], transVel_correction[2],
     angVel_correction[0], angVel_correction[1], angVel_correction[2]);
  #endif

  #pragma omp parallel for schedule(dynamic)
  for(size_t i=0; i<vInfo.size(); i++) {
    BlockInfo info = vInfo[i];
    const auto pos = obstacleBlocks.find(info.blockID);
    if(pos == obstacleBlocks.end()) continue;

    for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
    for(int iy=0; iy<FluidBlock::sizeY; ++iy)
    for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
      double p[3];
      info.pos(p, ix, iy, iz);
      p[0]-=CoM[0];
      p[1]-=CoM[1];
      p[2]-=CoM[2];
      const double correctVel[3] = {
        transVel_correction[0] +(angVel_correction[1]*p[2] - angVel_correction[2]*p[1]),
        transVel_correction[1] +(angVel_correction[2]*p[0] - angVel_correction[0]*p[2]),
        transVel_correction[2] +(angVel_correction[0]*p[1] - angVel_correction[1]*p[0])
      };
      pos->second->udef[iz][iy][ix][0] -= correctVel[0];
      pos->second->udef[iz][iy][ix][1] -= correctVel[1];
      pos->second->udef[iz][iy][ix][2] -= correctVel[2];
    }
  }
  #ifndef NDEBUG
    double dummy_ang[3], dummy_lin[3];
    _computeUdefMoments(dummy_lin, dummy_ang, CoM);
    const double EPS = 10*std::numeric_limits<Real>::epsilon();
    assert(std::fabs(dummy_lin[0])<EPS);
    assert(std::fabs(dummy_lin[1])<EPS);
    assert(std::fabs(dummy_lin[1])<EPS);
    assert(std::fabs(dummy_ang[0])<EPS);
    assert(std::fabs(dummy_ang[1])<EPS);
    assert(std::fabs(dummy_ang[2])<EPS);
  #endif
}

void IF3D_ObstacleOperator::_parseArguments(ArgumentParser & parser)
{
    parser.set_strict_mode();
    length = parser("-L").asDouble();
    position[0] = parser("-xpos").asDouble();
    parser.unset_strict_mode();
    position[1] = parser("-ypos").asDouble(ext_Y/2);
    position[2] = parser("-zpos").asDouble(ext_Z/2);
    quaternion[0] = parser("-quat0").asDouble(1.0);
    quaternion[1] = parser("-quat1").asDouble(0.0);
    quaternion[2] = parser("-quat2").asDouble(0.0);
    quaternion[3] = parser("-quat3").asDouble(0.0);
    _2Dangle = 2*std::atan2(quaternion[3], quaternion[0]);
    if(!rank)
    printf("Obstacle L=%g, pos=[%g %g %g], q=[%g %g %g %g]\n",
      length,position[0],position[1],position[2],quaternion[0],quaternion[1],quaternion[2],quaternion[3]);

    const double one = sqrt(quaternion[0]*quaternion[0]
         +quaternion[1]*quaternion[1]
         +quaternion[2]*quaternion[2]
         +quaternion[3]*quaternion[3]);

    if(fabs(one-1.0) > 5*numeric_limits<Real>::epsilon()) {
      printf("Parsed quaternion length is not equal to one. It really ought to be.\n");
      fflush(0);
      abort();
    }
    if(length < 5*numeric_limits<Real>::epsilon()) {
      printf("Parsed length is equal to zero. It really ought not to be.\n");
      fflush(0);
      abort();
    }

    // if true, obstacle will never change its velocity:
    // bForcedInLabFrame = parser("-bForcedInLabFrame").asBool(false);
    bool bFSM_alldir = parser("-bForcedInSimFrame").asBool(false);
    bForcedInSimFrame[0] = bFSM_alldir || parser("-bForcedInSimFrame_x").asBool(false);
    bForcedInSimFrame[1] = bFSM_alldir || parser("-bForcedInSimFrame_y").asBool(false);
    bForcedInSimFrame[2] = bFSM_alldir || parser("-bForcedInSimFrame_z").asBool(false);

    if(bForcedInSimFrame[0]) {
      const double xvel = -parser("-xvel").asDouble(0);
      transVel_imposed[0] = xvel; transVel[0] = xvel;
      if(!rank)
      printf("Obstacle forced to move relative to sim domain with constant x-vel:%f ", xvel);
    }
    if(bForcedInSimFrame[1]) {
      const double yvel = -parser("-yvel").asDouble(0);
      transVel_imposed[1] = yvel; transVel[1] = yvel;
      if(!rank)
      printf("Obstacle forced to move relative to sim domain with constant y-vel:%f ", yvel);
    }
    if(bForcedInSimFrame[2]) {
      const double zvel = -parser("-zvel").asDouble(0);
      transVel_imposed[2] = zvel; transVel[2] = zvel;
      if(!rank)
      printf("Obstacle forced to move relative to sim domain with constant z-vel:%f ", zvel);
    }
    bFixToPlanar = parser("-bFixToPlanar").asBool(false);


    const bool anyVelForced = bForcedInSimFrame[0] || bForcedInSimFrame[1] || bForcedInSimFrame[2];
    if(anyVelForced) {
      if(!rank) printf("and no angular velocity.\n");
      bBlockRotation[0] = true;
      bBlockRotation[1] = true;
      bBlockRotation[2] = true;
    }
    // this is different, obstacle can change the velocity, but sim frame will follow:
    bool bFOR_alldir = parser("-bFixFrameOfRef").asBool(false);
    bFixFrameOfRef[0] = bFOR_alldir || parser("-bFixFrameOfRef_x").asBool(false);
    bFixFrameOfRef[1] = bFOR_alldir || parser("-bFixFrameOfRef_y").asBool(false);
    bFixFrameOfRef[2] = bFOR_alldir || parser("-bFixFrameOfRef_z").asBool(false);
    bForces = parser("-computeForces").asBool(true);
}

void IF3D_ObstacleOperator::_writeComputedVelToFile(const int step_id, const double t, const Real* Uinf)
{
 if(rank!=0) return;
 stringstream ssR;
 ssR<<"computedVelocity_"<<obstacleID<<".dat";
    std::stringstream &savestream = logger.get_stream(ssR.str());
    const std::string tab("\t");

    if(step_id==0 && not printedHeaderVels) {
      printedHeaderVels = true;
      savestream<<"step"<<tab<<"time"<<tab<<"CMx"<<tab<<"CMy"<<tab<<"CMz"<<tab
      <<"quat_0"<<tab<<"quat_1"<<tab<<"quat_2"<<tab<<"quat_3"<<tab
      <<"vel_x"<<tab<<"vel_y"<<tab<<"vel_z"<<tab
      <<"angvel_x"<<tab<<"angvel_y"<<tab<<"angvel_z"<<tab<<"volume"<<tab
      <<"J0"<<tab<<"J1"<<tab<<"J2"<<tab<<"J3"<<tab<<"J4"<<tab<<"J5"<<std::endl;
    }

    savestream<<step_id<<tab;
    savestream.setf(std::ios::scientific);
    savestream.precision(std::numeric_limits<float>::digits10 + 1);
    savestream<<t<<tab<<position[0]<<tab<<position[1]<<tab<<position[2]<<tab
        <<quaternion[0]<<tab<<quaternion[1]<<tab<<quaternion[2]<<tab<<quaternion[3]<<tab
     <<transVel[0]-Uinf[0]<<tab<<transVel[1]-Uinf[1]<<tab<<transVel[2]-Uinf[2]<<tab
     <<angVel[0]<<tab<<angVel[1]<<tab<<angVel[2]<<tab<<volume<<tab
     <<J[0]<<tab<<J[1]<<tab<<J[2]<<tab<<J[3]<<tab<<J[4]<<tab<<J[5]<<std::endl;
}

void IF3D_ObstacleOperator::_writeDiagForcesToFile(const int step_id, const double t)
{
  if(rank!=0) return;
  stringstream ssR;
  ssR<<"diagnosticsForces_"<<obstacleID<<".dat";
  std::stringstream &savestream = logger.get_stream(ssR.str());
  const std::string tab("\t");

  if(step_id==0)
    savestream << "step" << tab << "time" << tab << "mass" << tab
     << "force_x" << tab << "force_y" << tab << "force_z" << tab
     << "torque_x" << tab << "torque_y" << tab << "torque_z" << std::endl;

  savestream << step_id << tab;
  savestream.setf(std::ios::scientific);
  savestream.precision(std::numeric_limits<float>::digits10 + 1);
  savestream<<t<<tab<<mass<<tab<<force[0]<<tab<<force[1]<<tab<<force[2]<<tab
      <<torque[0]<<tab<<torque[1]<<tab<<torque[2]<<std::endl;
}

void IF3D_ObstacleOperator::computeDiagnostics(const int stepID, const double time, const Real* Uinf, const double lambda)
{
  double CM[3];
  this->getCenterOfMass(CM);
  double _area=0,_forcex=0,_forcey=0,_forcez=0,_torquex=0,_torquey=0,_torquez=0;

  #pragma omp parallel for schedule(dynamic) reduction(+:_area,_forcex,_forcey,_forcez,_torquex,_torquey,_torquez)
  for(size_t i=0; i<vInfo.size(); i++) {
      BlockInfo info = vInfo[i];
      const auto pos = obstacleBlocks.find(info.blockID);
      if(pos == obstacleBlocks.end()) continue;
      FluidBlock& b = *(FluidBlock*)info.ptrBlock;

      for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
      for(int iy=0; iy<FluidBlock::sizeY; ++iy)
      for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
          const double Xs = pos->second->chi[iz][iy][ix];
          if (Xs <= 0) continue;
          double p[3];
          info.pos(p, ix, iy, iz);
          p[0]-=CM[0];
          p[1]-=CM[1];
          p[2]-=CM[2];
        const double object_UR[3] = {
          angVel[1]*p[2]-angVel[2]*p[1],
          angVel[2]*p[0]-angVel[0]*p[2],
          angVel[0]*p[1]-angVel[1]*p[0]
        };
          const double object_UDEF[3] = {
              pos->second->udef[iz][iy][ix][0],
              pos->second->udef[iz][iy][ix][1],
              pos->second->udef[iz][iy][ix][2]
          };
          const double U[3] = {
              b(ix,iy,iz).u+Uinf[0] - (transVel[0]+object_UR[0]+object_UDEF[0]),
              b(ix,iy,iz).v+Uinf[1] - (transVel[1]+object_UR[1]+object_UDEF[1]),
              b(ix,iy,iz).w+Uinf[2] - (transVel[2]+object_UR[2]+object_UDEF[2])
          };
          _area += Xs;
          _forcex += U[0]*Xs;
          _forcey += U[1]*Xs;
          _forcez += U[2]*Xs;
          _torquex += (p[1]*U[2]-p[2]*U[1])*Xs;
          _torquey += (p[2]*U[0]-p[0]*U[2])*Xs;
          _torquez += (p[0]*U[1]-p[1]*U[0])*Xs;
      }
  }

  double globals[7] = {0,0,0,0,0,0,0};
  double locals[7] = {_forcex,_forcey,_forcez,_torquex,_torquey,_torquez,_area};
  MPI_Allreduce(locals, globals, 7, MPI_DOUBLE, MPI_SUM, grid->getCartComm());

  const double dV = std::pow(vInfo[0].h_gridpoint, 3);
  force[0]  = globals[0]*dV*lambda;
  force[1]  = globals[1]*dV*lambda;
  force[2]  = globals[2]*dV*lambda;
  torque[0] = globals[3]*dV*lambda;
  torque[1] = globals[4]*dV*lambda;
  torque[2] = globals[5]*dV*lambda;
  mass      = globals[6]*dV;
  #ifndef __RL_TRAINING
  _writeDiagForcesToFile(stepID, time);
  #endif
}

void IF3D_ObstacleOperator::computeVelocities(const Real* Uinf)
{
  double CM[3];
  this->getCenterOfMass(CM);
  const double h  = vInfo[0].h_gridpoint;
  const double dv = std::pow(vInfo[0].h_gridpoint,3);
  {
    double V = 0, lm0 = 0, lm1 = 0, lm2 = 0; //linear momenta
    #pragma omp parallel for schedule(dynamic) reduction(+:V,lm0,lm1,lm2)
    for(size_t i=0; i<vInfo.size(); i++) {
      BlockInfo info = vInfo[i];
      FluidBlock & b = *(FluidBlock*)info.ptrBlock;
      const auto pos = obstacleBlocks.find(info.blockID);
      if(pos == obstacleBlocks.end()) continue;

      for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
      for(int iy=0; iy<FluidBlock::sizeY; ++iy)
      for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
        const double Xs = pos->second->chi[iz][iy][ix];
        if (Xs <= 0) continue;
        V     += Xs;
        lm0   += Xs * b(ix,iy,iz).u;
        lm1   += Xs * b(ix,iy,iz).v;
        lm2   += Xs * b(ix,iy,iz).w;
      }
    }

    double globals[4] = {0,0,0,0};
    double locals[4] = {lm0,lm1,lm2,V};
    MPI_Allreduce(locals, globals, 4, MPI_DOUBLE, MPI_SUM, grid->getCartComm());
    assert(globals[3] > std::numeric_limits<double>::epsilon());

    if(bForcedInSimFrame[0]) {
      transVel[0] = transVel_imposed[0];
      transVel_computed[0] = globals[0]/globals[3] + Uinf[0];
    } else
      transVel[0]          = globals[0]/globals[3] + Uinf[0];
    if(bForcedInSimFrame[1]) {
      transVel[1] = transVel_imposed[1];
      transVel_computed[1] = globals[1]/globals[3] + Uinf[1];
    } else
      transVel[1]          = globals[1]/globals[3] + Uinf[1];
    if(bForcedInSimFrame[2]) {
      transVel[2] = transVel_imposed[2];
      transVel_computed[2] = globals[2]/globals[3] + Uinf[2];
    } else
      transVel[2]          = globals[2]/globals[3] + Uinf[2];
    volume      = globals[3] * dv;
    if(bFixToPlanar && not bForcedInSimFrame[2]) transVel[2] = 0.0;
  }
  {
    double J0=0, J1=0, J2=0, J3=0, J4=0, J5=0, am0=0, am1=0, am2=0;
    #pragma omp parallel for schedule(dynamic) reduction(+:J0,J1,J2,J3,J4,J5,am0,am1,am2)
    for(size_t i=0; i<vInfo.size(); i++) {
      BlockInfo info = vInfo[i];
      FluidBlock & b = *(FluidBlock*)info.ptrBlock;
      const auto pos = obstacleBlocks.find(info.blockID);
      if(pos == obstacleBlocks.end()) continue;

      for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
      for(int iy=0; iy<FluidBlock::sizeY; ++iy)
      for(int ix=0; ix<FluidBlock::sizeX; ++ix) {
        const double Xs = pos->second->chi[iz][iy][ix];
        if (Xs <= 0) continue;
        double p[3];
        info.pos(p, ix, iy, iz);
        p[0]-=CM[0];
        p[1]-=CM[1];
        p[2]-=CM[2];
        const double u_= b(ix,iy,iz).u, v_= b(ix,iy,iz).v, w_= b(ix,iy,iz).w;

        am0 += Xs * (p[1]*(w_-transVel[2]) - p[2]*(v_-transVel[1]));
        am1 += Xs * (p[2]*(u_-transVel[0]) - p[0]*(w_-transVel[2]));
        am2 += Xs * (p[0]*(v_-transVel[1]) - p[1]*(u_-transVel[0]));
        J0  += Xs * (p[1]*p[1]+p[2]*p[2]);
        J1  += Xs * (p[0]*p[0]+p[2]*p[2]);
        J2  += Xs * (p[0]*p[0]+p[1]*p[1]);
        J3  -= Xs *  p[0]*p[1];
        J4  -= Xs *  p[0]*p[2];
        J5  -= Xs *  p[1]*p[2];
      }
    }

    double globals[9] = {0,0,0,0,0,0,0,0,0};
    double locals[9] = {am0,am1,am2,J0,J1,J2,J3,J4,J5};
    MPI_Allreduce(locals, globals, 9, MPI_DOUBLE, MPI_SUM, grid->getCartComm());

    if(bFixToPlanar) {
      angVel[0] = angVel[1] = 0.0;
      (bBlockRotation[2]? angVel_computed[2] : angVel[2]) = globals[2]/globals[5]; // av2/j2
    } else {
      //solve avel = invJ \dot angMomentum, do not multiply by h^3 for numerics
      const double AM[3] = {globals[0]*h, globals[1]*h, globals[2]*h};
      const double J_[6] = {globals[3]*h, globals[4]*h, globals[5]*h,
                          globals[6]*h, globals[7]*h, globals[8]*h};
      const double detJ = J_[0]*(J_[1]*J_[2] - J_[5]*J_[5])+
                        J_[3]*(J_[4]*J_[5] - J_[2]*J_[3])+
                        J_[4]*(J_[3]*J_[5] - J_[1]*J_[4]);
      const double invDetJ = 1./detJ;
      const double invJ[6] = {
        invDetJ * (J_[1]*J_[2] - J_[5]*J_[5]),
        invDetJ * (J_[0]*J_[2] - J_[4]*J_[4]),
        invDetJ * (J_[0]*J_[1] - J_[3]*J_[3]),
        invDetJ * (J_[4]*J_[5] - J_[2]*J_[3]),
        invDetJ * (J_[3]*J_[5] - J_[1]*J_[4]),
        invDetJ * (J_[3]*J_[4] - J_[0]*J_[5])
      };
      (bBlockRotation[0]? angVel_computed[0] : angVel[0]) = invJ[0]*AM[0] + invJ[3]*AM[1] + invJ[4]*AM[2];
      (bBlockRotation[1]? angVel_computed[1] : angVel[1]) = invJ[3]*AM[0] + invJ[1]*AM[1] + invJ[5]*AM[2];
      (bBlockRotation[2]? angVel_computed[2] : angVel[2]) = invJ[4]*AM[0] + invJ[5]*AM[1] + invJ[2]*AM[2];
    }
    J[0] = globals[3] * dv;
    J[1] = globals[4] * dv;
    J[2] = globals[5] * dv;
    J[3] = globals[6] * dv;
    J[4] = globals[7] * dv;
    J[5] = globals[8] * dv;
  }
}

void IF3D_ObstacleOperator::dumpWake(const int stepID, const double t, const Real* Uinf)
{
  //horrible, dont look at it!!!
  stringstream ssR;
  ssR<<"wakeValues_rank"<<rank<<"ID"<<obstacleID<<"_"<<t<<".dat";
  FILE * pFile = fopen (ssR.str().c_str(), "ab");
  DumpWake kernel(Uinf, position, pFile, length);
  SynchronizerMPI& Synch = grid->sync(kernel);
  LabMPI labs;
  labs.prepare(*grid, Synch);
  MPI_Barrier(grid->getCartComm());
  vector<BlockInfo> avail0 = Synch.avail_inner();
  vector<BlockInfo> avail1 = Synch.avail_halo();
  std::map<int, BlockInfo*> tmp;
  for(size_t i=0; i<vInfo.size(); i++)
  {
    const int look4 = vInfo[i].blockID;
    bool found = false;
    for (size_t j=0; j<avail0.size(); j++) {
    if(avail0[j].blockID == look4) {
        if(found) printf("Two blocks with the same ID?!\n");
        else tmp[look4] = &avail0[j];
        found = true;
      }
    }
    for (size_t j=0; j<avail1.size(); j++) {
      if(avail1[j].blockID == look4) {
        if(found) printf("Two blocks with the same ID?!\n");
        else tmp[look4] = &avail1[j];
        found = true;
      }
    }
    if(!found) printf("Wtf missing blocks?? Brace for segfault!\n");
  }
  //now in tmp i have addresses to all info, and i can dump with same sorting as vInfo
  for(size_t i=0; i<vInfo.size(); i++)
  { //i care more about ease of postprocess here than scaling
    const int blockID = vInfo[i].blockID;
    BlockInfo info = *tmp[blockID];
    FluidBlock& b = *(FluidBlock*)info.ptrBlock;
    labs.load(info, 0);
    kernel(labs, info, b);
  }
  MPI_Barrier(grid->getCartComm());
  fclose (pFile);
  /*
  if(!rank)
  {
    stringstream ssR;
    ssR<<"headValues_ID"<<obstacleID<<"_"<<t<<".dat";
    FILE * pFile = fopen (ssR.str().c_str(), "ab");
    assert(sr.VelNAbove.size()==sr.NpLatLine);
    for(int i=sr.NpLatLine-1; i>=0; --i) {
      const double d[] = {
        sr.PXAbove[i]-position[0],  sr.PYAbove[i]-position[1],
        sr.VelNAbove[i],            sr.VelTAbove[i]
      };
      fwrite(d,sizeof(double),4,pFile);
    }
    for(int i=0; i<sr.NpLatLine; ++i){
      const double d[] = {
        sr.PXBelow[i]-position[0],  sr.PYBelow[i]-position[1],
        sr.VelNBelow[i],            sr.VelTBelow[i]
      };
      fwrite(d,sizeof(double),4,pFile);
    }
    fclose (pFile);
  }
  */
}

void IF3D_ObstacleOperator::computeForces(const int stepID, const double time,
  const double dt, const Real* Uinf, const double NU, const bool bDump)
{
  //TODO: improve dumping: gather arrays before writing to file
  if(not bForces) return;

  double CM[3];
  this->getCenterOfMass(CM);
  const double velx_tot = transVel[0]-Uinf[0];
  const double vely_tot = transVel[1]-Uinf[1];
  const double velz_tot = transVel[2]-Uinf[2];

  #ifdef __RL_TRAINING
    if(!bInteractive) {
      sr.updateAverages(dt,_2Dangle,velx_tot,vely_tot,angVel[2],0,0,0,0,0,0,0,0,0,0);
      return;
    }
  #endif

  double vel_unit[3] = {0., 0., 0.};
  const double vel_norm = std::sqrt(velx_tot*velx_tot
                                + vely_tot*vely_tot
                                + velz_tot*velz_tot);
  if (vel_norm>1e-9) {
      vel_unit[0] = velx_tot/vel_norm;
      vel_unit[1] = vely_tot/vel_norm;
      vel_unit[2] = velz_tot/vel_norm;
  }
  const int nthreads = omp_get_max_threads();
  vector<OperatorComputeForces*> forces(nthreads, nullptr);
  for(int i=0;i<nthreads;++i)
    forces[i] = new OperatorComputeForces(NU,dt,vel_unit,Uinf,CM);

  compute<OperatorComputeForces, SURFACE>(forces);

  for(int i=0; i<nthreads; i++) delete forces[i];

  static const int nQoI = ObstacleBlock::nQoI;
  sum = vector<double>(nQoI, 0);
  for (auto & block : obstacleBlocks) block.second->sumQoI(sum);

  MPI_Allreduce(MPI_IN_PLACE, sum.data(), nQoI, MPI_DOUBLE, MPI_SUM, grid->getCartComm());

  //additive quantities: (check against order in sumQoI of ObstacleBlocks.h )
  unsigned k = 0;
  surfForce[0]= sum[k++]; surfForce[1]= sum[k++]; surfForce[2]= sum[k++];
  double forcex_P= sum[k++], forcey_P= sum[k++], forcez_P= sum[k++];
  double forcex_V= sum[k++], forcey_V= sum[k++], forcez_V= sum[k++];
  double torquex= sum[k++], torquey= sum[k++], torquez= sum[k++];
  double gammax= sum[k++], gammay= sum[k++], gammaz= sum[k++];

  drag = sum[k++]; thrust = sum[k++]; Pout = sum[k++];
  PoutBnd = sum[k++]; defPower = sum[k++]; defPowerBnd = sum[k++];

  //derived quantities:
  Pthrust    = thrust*vel_norm;
  Pdrag      =   drag*vel_norm;
  const double eps = std::numeric_limits<Real>::epsilon();
  EffPDef    = Pthrust/(Pthrust-min(defPower,(double)0.)+eps);
  EffPDefBnd = Pthrust/(Pthrust-    defPowerBnd         +eps);

  sr.updateAverages(dt,_2Dangle, velx_tot, vely_tot, angVel[2], Pout, PoutBnd,
    defPower, defPowerBnd, EffPDef, EffPDefBnd, Pthrust, Pdrag, thrust, drag);

  #ifndef __RL_TRAINING
  if (bDump) {
    char buf[500];
    sprintf(buf, "surface_%02d_%07d_rank%03d.raw", obstacleID, stepID, rank);
    FILE * pFile = fopen (buf, "wb");
    for(auto & block : obstacleBlocks) block.second->print(pFile);
    fflush(pFile);
    fclose(pFile);
  }

  if(rank==0) {
    stringstream ssF, ssP;
    ssF<<"forceValues_"<<obstacleID<<".dat";
    ssP<<"powerValues_"<<obstacleID<<".dat";

    std::stringstream &fileForce = logger.get_stream(ssF.str());
    if(stepID==0)
      fileForce<<"time Fx Fy Fz FxPres FyPres FzPres FxVisc FyVisc FzVisc TorqX TorqY TorqZ Gx Gy Gz drag thrust"<<std::endl;

    fileForce<<time<<" "<<surfForce[0]<<" "<<surfForce[1]<<" "<<surfForce[2]
      <<" "<<forcex_P<<" "<<forcey_P<<" "<<forcez_P<<" "<<forcex_V<<" "
      <<forcey_V<<" "<<forcez_V<<" "<<torquex<<" "<<torquey<<" "<<torquez
      <<" "<<gammax<<" "<< gammay<<" "<<gammaz<<" "<<drag<<" "<<thrust<<endl;

    std::stringstream &filePower = logger.get_stream(ssP.str());
    if(stepID==0)
      filePower<<"time Pthrust Pdrag PoutBnd Pout defPowerBnd defPower EffPDefBnd EffPDef"<<std::endl;

    filePower<<time<<" "<<Pthrust<<" "<<Pdrag<<" "<<PoutBnd<<" "<<Pout<<" "
      <<defPowerBnd<<" "<<defPower<<" "<<EffPDefBnd<<" "<<EffPDef<<endl;
  }
  #endif
}

void IF3D_ObstacleOperator::update(const int step_id, const double t, const double dt, const Real* Uinf)
{
  position[0] += dt*transVel[0];
  position[1] += dt*transVel[1];
  position[2] += dt*transVel[2];

  const double dqdt[4] = {
    .5*( -angVel[0]*quaternion[1]-angVel[1]*quaternion[2]-angVel[2]*quaternion[3]),
    .5*(  angVel[0]*quaternion[0]+angVel[1]*quaternion[3]-angVel[2]*quaternion[2]),
    .5*( -angVel[0]*quaternion[3]+angVel[1]*quaternion[0]+angVel[2]*quaternion[1]),
    .5*(  angVel[0]*quaternion[2]-angVel[1]*quaternion[1]+angVel[2]*quaternion[0])
  };

  // normality preserving advection (Simulation of colliding constrained rigid bodies - Kleppmann 2007 Cambridge University, p51)
  // move the correct distance on the quaternion unit ball surface, end up with normalized quaternion
  const double deltaq[4] = {
      dqdt[0]*dt,
      dqdt[1]*dt,
      dqdt[2]*dt,
      dqdt[3]*dt
  };

  const double deltaq_length = std::sqrt(deltaq[0]*deltaq[0]+deltaq[1]*deltaq[1]+deltaq[2]*deltaq[2]+deltaq[3]*deltaq[3]);

  if(deltaq_length>std::numeric_limits<double>::epsilon()) {
      const double tanfac = std::tan(deltaq_length)/deltaq_length;
      const double num[4] = {
          quaternion[0]+tanfac*deltaq[0],
          quaternion[1]+tanfac*deltaq[1],
          quaternion[2]+tanfac*deltaq[2],
          quaternion[3]+tanfac*deltaq[3],
      };

      const double invDenum = 1./(std::sqrt(num[0]*num[0]+num[1]*num[1]+num[2]*num[2]+num[3]*num[3]));
      quaternion[0] = num[0]*invDenum;
      quaternion[1] = num[1]*invDenum;
      quaternion[2] = num[2]*invDenum;
      quaternion[3] = num[3]*invDenum;
  }

  //_2Dangle += dt*angVel[2];
  const double old2DA = _2Dangle;
  //keep consistency: get 2d angle from quaternions:
  _2Dangle = 2*std::atan2(quaternion[3], quaternion[0]);
  const double err = std::fabs(_2Dangle-old2DA-dt*angVel[2]);
  if(err>std::numeric_limits<Real>::epsilon() && !rank)
    printf("Discrepancy in angvel from quaternions: %f (%f %f)\n",
      err, (_2Dangle-old2DA)/dt, angVel[2]);

  double velx_tot = transVel[0]-Uinf[0];
  double vely_tot = transVel[1]-Uinf[1];
  double velz_tot = transVel[2]-Uinf[2];
  absPos[0] += dt*velx_tot;
  absPos[1] += dt*vely_tot;
  absPos[2] += dt*velz_tot;
  sr.updateInstant(position[0], absPos[0], position[1], absPos[1],
                    _2Dangle, velx_tot, vely_tot, angVel[2]);

  #ifndef NDEBUG
  if(rank==0) {
    #ifdef _VERBOSE_
     std::cout<<"POSITION INFO AFTER UPDATE T, DT: "<<t<<" "<<dt<<std::endl;
     std::cout<<"POS: "<<position[0]<<" "<<position[1]<<" "<<position[2]<<std::endl;
     std::cout<<"TVL: "<<transVel[0]<<" "<<transVel[1]<<" "<<transVel[2]<<std::endl;
     std::cout<<"QUT: "<<quaternion[0]<<" "<<quaternion[1]<<" "<<quaternion[2]<<" "<<quaternion[3]<<std::endl;
     std::cout<<"AVL: "<<angVel[0]<<" "<<angVel[1]<<" "<<angVel[2]<<std::endl;
    #endif
  }
  const double q_length=std::sqrt(quaternion[0]*quaternion[0]
        +  quaternion[1]*quaternion[1]
        +  quaternion[2]*quaternion[2]
        +  quaternion[3]*quaternion[3]);
  assert(std::abs(q_length-1.0) < 5*std::numeric_limits<Real>::epsilon());
  #endif
  #ifndef __RL_TRAINING
  _writeComputedVelToFile(step_id, t, Uinf);
  #endif
}

void IF3D_ObstacleOperator::characteristic_function()
{
  #pragma omp parallel
 {
  #pragma omp for schedule(dynamic)
  for(size_t i=0; i<vInfo.size(); i++) {
   BlockInfo info = vInfo[i];
   std::map<int,ObstacleBlock* >::const_iterator pos = obstacleBlocks.find(info.blockID);

   if(pos != obstacleBlocks.end()) {
        FluidBlock& b = *(FluidBlock*)info.ptrBlock;
    for(int iz=0; iz<FluidBlock::sizeZ; iz++)
    for(int iy=0; iy<FluidBlock::sizeY; iy++)
    for(int ix=0; ix<FluidBlock::sizeX; ix++)
     b(ix,iy,iz).chi = std::max(pos->second->chi[iz][iy][ix], b(ix,iy,iz).chi);
   }
  }
 }
}

std::vector<int> IF3D_ObstacleOperator::intersectingBlockIDs(const int buffer) const
{
 assert(buffer <= 2); // only works for 2: if different definition of deformation blocks, implement your own
 std::vector<int> retval;
 const int N = vInfo.size();

 for(int i=0; i<N; i++) {
  BlockInfo info = vInfo[i];
  std::map<int,ObstacleBlock* >::const_iterator pos = obstacleBlocks.find(info.blockID);
  if(pos != obstacleBlocks.end()) retval.push_back(info.blockID);
 }
 return retval;
}

void IF3D_ObstacleOperator::getSkinsAndPOV(Real& x, Real& y, Real& th,
  Real*& pXL, Real*& pYL, Real*& pXU, Real*& pYU, int& Npts)
{
  printf("Entered the wrong get skin operator\n");
 fflush(0);
  abort();
}

void IF3D_ObstacleOperator::interpolateOnSkin(const double time, const int stepID, bool dumpWake)
{
  //printf("Entered the wrong interpolate operator\n");
 //fflush(0);
  //abort();
}

void IF3D_ObstacleOperator::execute(const int iAgent, const double time, const vector<double> action)
{
  printf("Entered the wrong execute operator\n");
 fflush(0);
  abort();
}

void IF3D_ObstacleOperator::create(const int step_id,const double time,
    const double dt, const Real *Uinf)
{
  printf("Entered the wrong create operator\n");
 fflush(0);
  abort();
}

void IF3D_ObstacleOperator::computeChi(const int step_id, const double time,
  const double dt, const Real *Uinf, int& mpi_status) {}

void IF3D_ObstacleOperator::finalize(const int step_id,const double time,
  const double dt, const Real *Uinf)
{ }

void IF3D_ObstacleOperator::getTranslationVelocity(double UT[3]) const
{
    UT[0]=transVel[0];
    UT[1]=transVel[1];
    UT[2]=transVel[2];
}

void IF3D_ObstacleOperator::setTranslationVelocity(double UT[3])
{
    transVel[0] = UT[0];
    transVel[1] = UT[1];
    transVel[2] = UT[2];
}

void IF3D_ObstacleOperator::getAngularVelocity(double W[3]) const
{
    W[0]=angVel[0];
    W[1]=angVel[1];
    W[2]=angVel[2];
}

void IF3D_ObstacleOperator::setAngularVelocity(const double W[3])
{
   angVel[0]=W[0];
   angVel[1]=W[1];
   angVel[2]=W[2];
}

void IF3D_ObstacleOperator::getCenterOfMass(double CM[3]) const
{
    CM[0]=position[0];
    CM[1]=position[1];
    CM[2]=position[2];
}

void IF3D_ObstacleOperator::save(const int step_id, const double t, std::string filename)
{
    if(rank!=0) return;
    sr.save(step_id,filename);
    std::ofstream savestream;
    savestream.setf(std::ios::scientific);
    savestream.precision(std::numeric_limits<Real>::digits10 + 1);
    savestream.open(filename+".txt");
    savestream<<t<<std::endl;
    savestream<<position[0]<<"\t"<<position[1]<<"\t"<<position[2]<<std::endl;
    savestream<<absPos[0]<<"\t"<<absPos[1]<<"\t"<<absPos[2]<<std::endl;
    savestream<<quaternion[0]<<"\t"<<quaternion[1]<<"\t"<<quaternion[2]<<"\t"<<quaternion[3]<<std::endl;
    savestream<<transVel[0]<<"\t"<<transVel[1]<<"\t"<<transVel[2]<<std::endl;
    savestream<<angVel[0]<<"\t"<<angVel[1]<<"\t"<<angVel[2]<<std::endl;
    savestream<<_2Dangle<<std::endl;
}

void IF3D_ObstacleOperator::restart(const double t, std::string filename)
{
    sr.restart(filename);
    std::ifstream restartstream;
    restartstream.open(filename+".txt");
    if(!restartstream.good()){
      printf("Could not restart from file\n");
      return;
    }
    Real restart_time;
    restartstream >> restart_time;
    //assert(std::abs(restart_time-t) < 1e-9);

    restartstream>>position[0]>>position[1]>>position[2];
    restartstream>>absPos[0]>>absPos[1]>>absPos[2];
    restartstream>>quaternion[0]>>quaternion[1]>>quaternion[2]>>quaternion[3];
    restartstream>>transVel[0]>>transVel[1]>>transVel[2];
    restartstream>>angVel[0]>>angVel[1]>>angVel[2];
    restartstream >> _2Dangle;
    restartstream.close();

    {
        std::cout<<"RESTARTED BODY: "<<std::endl;
        std::cout<<"TIME: \t"<<restart_time<<std::endl;
        std::cout<<"POS : \t"<<position[0]<<" "<<position[1]
                               <<" "<<position[2]<<std::endl;
        std::cout<<"ABS POS : \t"<<absPos[0]<<" "<<absPos[1]
                               <<" "<<absPos[2]<<std::endl;
        std::cout<<"ANGLE:\t"<<quaternion[0]<<" "<<quaternion[1]<<" "
                               <<quaternion[2]<<" "<<quaternion[3]<<std::endl;
        std::cout<<"TVEL: \t"<<transVel[0]<<" "<<transVel[1]
                               <<" "<<transVel[2]<<std::endl;
        std::cout<<"AVEL: \t"<<angVel[0]<<" "<<angVel[1]
                               <<" "<<angVel[2]<<std::endl;
        std::cout<<"2D angle: \t"<<_2Dangle<<std::endl;
    }
}

void IF3D_ObstacleOperator::Accept(ObstacleVisitor * visitor)
{
 visitor->visit(this);
}
