//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#include "Obstacle.h"
#include "../Utils/BufferedLogger.h"

#include <Cubism/ArgumentParser.h>
#include <gsl/gsl_linalg.h>
#include <fstream>

CubismUP_3D_NAMESPACE_BEGIN
using namespace cubism;

using UDEFMAT = Real[CUP_BLOCK_SIZEZ][CUP_BLOCK_SIZEY][CUP_BLOCK_SIZEX][3];
using CHIMAT =  Real[CUP_BLOCK_SIZEZ][CUP_BLOCK_SIZEY][CUP_BLOCK_SIZEX];
static constexpr Real EPS = std::numeric_limits<Real>::epsilon();

Obstacle::Obstacle(SimulationData& s,  ArgumentParser& parser): sim(s)
{
  length = parser("-L").asDouble();          // Mandatory.
  position[0] = parser("-xpos").asDouble();  // Mandatory.
  position[1] = parser("-ypos").asDouble(sim.extents[1] / 2);
  position[2] = parser("-zpos").asDouble(sim.extents[2] / 2);
  quaternion[0] = parser("-quat0").asDouble(0.0);
  quaternion[1] = parser("-quat1").asDouble(0.0);
  quaternion[2] = parser("-quat2").asDouble(0.0);
  quaternion[3] = parser("-quat3").asDouble(0.0);

  Real planarAngle = parser("-planarAngle").asDouble(0.0) / 180 * M_PI;
  const Real q_length = std::sqrt(quaternion[0]*quaternion[0] + quaternion[1]*quaternion[1] 
                                + quaternion[2]*quaternion[2] + quaternion[3]*quaternion[3]);
  quaternion[0] /= q_length;
  quaternion[1] /= q_length;
  quaternion[2] /= q_length;
  quaternion[3] /= q_length;
  if(std::fabs(q_length-1.0) > 100*EPS)
  {
    quaternion[0] = std::cos(0.5*planarAngle);
    quaternion[1] = 0;
    quaternion[2] = 0;
    quaternion[3] = std::sin(0.5*planarAngle);
  } 
  else
  {
    if(std::fabs(planarAngle) > 0 && sim.rank == 0)
      std::cout << "WARNING: Obstacle arguments include both quaternions and planarAngle." <<
                   "Quaterion arguments have priority and therefore planarAngle will be ignored.\n";
    planarAngle = 2 * std::atan2(quaternion[3], quaternion[0]);
  }

  // if true, obstacle will never change its velocity:
  bool bFSM_alldir = parser("-bForcedInSimFrame").asBool(false);
  bForcedInSimFrame[0] = bFSM_alldir || parser("-bForcedInSimFrame_x").asBool(false);
  bForcedInSimFrame[1] = bFSM_alldir || parser("-bForcedInSimFrame_y").asBool(false);
  bForcedInSimFrame[2] = bFSM_alldir || parser("-bForcedInSimFrame_z").asBool(false);
  // only active if corresponding bForcedInLabFrame is true:
  Real enforcedVelocity [3];
  enforcedVelocity[0] = -parser("-xvel").asDouble(0.0);
  enforcedVelocity[1] = -parser("-yvel").asDouble(0.0);
  enforcedVelocity[2] = -parser("-zvel").asDouble(0.0);
  const bool bFixToPlanar = parser("-bFixToPlanar").asBool(false);
  // this is different, obstacle can change the velocity, but sim frame will follow:
  bool bFOR_alldir = parser("-bFixFrameOfRef").asBool(false);
  bFixFrameOfRef[0] = bFOR_alldir || parser("-bFixFrameOfRef_x").asBool(false);
  bFixFrameOfRef[1] = bFOR_alldir || parser("-bFixFrameOfRef_y").asBool(false);
  bFixFrameOfRef[2] = bFOR_alldir || parser("-bFixFrameOfRef_z").asBool(false);
  // boolean to break symmetry to trigger vortex shedding
  bBreakSymmetry = parser("-bBreakSymmetry").asBool(false);

  absPos[0] = position[0]; 
  absPos[1] = position[1]; 
  absPos[2] = position[2];

  if (!sim.rank)
  {
    printf("Obstacle L=%g, pos=[%g %g %g], q=[%g %g %g %g]\n",
           length, position[0], position[1], position[2],
           quaternion[0], quaternion[1], quaternion[2], quaternion[3]);
  }

  const Real one = std::sqrt(
          quaternion[0] * quaternion[0] + quaternion[1] * quaternion[1]
        + quaternion[2] * quaternion[2] + quaternion[3] * quaternion[3]);

  if (std::fabs(one - 1.0) > 5 * EPS)
  {
    printf("Parsed quaternion length is not equal to one. It really ought to be.\n");
    fflush(0);
    abort();
  }
  if (length < 5 * EPS)
  {
    printf("Parsed length is equal to zero. It really ought not to be.\n");
    fflush(0);
    abort();
  }

  for (int d = 0; d < 3; ++d)
  {
    bForcedInSimFrame[d] = bForcedInSimFrame[d];
    if (bForcedInSimFrame[d])
    {
      transVel_imposed[d] = transVel[d] = enforcedVelocity[d];
      if (!sim.rank)
        printf("Obstacle forced to move relative to sim domain with constant %c-vel: %f\n","xyz"[d], transVel[d]);
    }
  }

  const bool anyVelForced = bForcedInSimFrame[0] || bForcedInSimFrame[1] || bForcedInSimFrame[2];
  if(anyVelForced)
  {
    if (!sim.rank) printf("Obstacle has no angular velocity.\n");
    bBlockRotation[0] = true;
    bBlockRotation[1] = true;
    bBlockRotation[2] = true;
  }

  if(bFixToPlanar)
  {
    if (!sim.rank) printf("Obstacle motion restricted to constant Z-plane.\n");
    bForcedInSimFrame[2] = true;
    transVel_imposed[2] = 0;
    bBlockRotation[1] = true;
    bBlockRotation[0] = true;
  }

  if( bBreakSymmetry )
    if (!sim.rank) printf("Symmetry broken by imposing sinusodial y-velocity in t=[1,2].\n");
}

void Obstacle::updateLabVelocity( int nSum[3], Real uSum[3] )
{
  if (bFixFrameOfRef[0]) { nSum[0]+=1; uSum[0] -= transVel[0]; }
  if (bFixFrameOfRef[1]) { nSum[1]+=1; uSum[1] -= transVel[1]; }
  if (bFixFrameOfRef[2]) { nSum[2]+=1; uSum[2] -= transVel[2]; }
}

void Obstacle::computeVelocities()
{
  std::vector<double> A(36);//need to use double (not Real) for GSL
  A[0*6 + 0] =      penalM ; A[0*6 + 1] =         0.0 ; A[0*6 + 2] =         0.0; A[0*6 + 3] =         0.0; A[0*6 + 4] = +penalCM[2]; A[0*6 + 5] = -penalCM[1];
  A[1*6 + 0] =         0.0 ; A[1*6 + 1] =      penalM ; A[1*6 + 2] =         0.0; A[1*6 + 3] = -penalCM[2]; A[1*6 + 4] =         0.0; A[1*6 + 5] = +penalCM[0];
  A[2*6 + 0] =         0.0 ; A[2*6 + 1] =         0.0 ; A[2*6 + 2] =      penalM; A[2*6 + 3] = +penalCM[1]; A[2*6 + 4] = -penalCM[0]; A[2*6 + 5] =         0.0;
  A[3*6 + 0] =         0.0 ; A[3*6 + 1] = -penalCM[2] ; A[3*6 + 2] = +penalCM[1]; A[3*6 + 3] =   penalJ[0]; A[3*6 + 4] =   penalJ[3]; A[3*6 + 5] =   penalJ[4];
  A[4*6 + 0] = +penalCM[2] ; A[4*6 + 1] =         0.0 ; A[4*6 + 2] = -penalCM[0]; A[4*6 + 3] =   penalJ[3]; A[4*6 + 4] =   penalJ[1]; A[4*6 + 5] =   penalJ[5];
  A[5*6 + 0] = -penalCM[1] ; A[5*6 + 1] = +penalCM[0] ; A[5*6 + 2] =         0.0; A[5*6 + 3] =   penalJ[4]; A[5*6 + 4] =   penalJ[5]; A[5*6 + 5] =   penalJ[2];

  // TODO here we can add dt * appliedForce/Torque[i]
  double b[6] = { //need to use double (not Real) for GSL
    penalLmom[0], penalLmom[1], penalLmom[2],
    penalAmom[0], penalAmom[1], penalAmom[2]
  };

  // modify y-velocity for symmetry breaking
  if( bBreakSymmetry )
  {
    if( sim.time>3.0 && sim.time<4.0 )
      // transVel_imposed[1] = length*std::sin(M_PI*(sim.time-3.0)); // for Re=300
      transVel_imposed[1] = 0.1*length*std::sin(M_PI*(sim.time-3.0)); // for Re=1000
    else
      transVel_imposed[1] = 0.0;
  }

  //Momenta are conserved if a dof (a row of mat A) is not externally forced
  //This means that if obstacle is free to move according to fluid forces,
  //momenta after penal should be equal to moments before penal!
  //If dof is forced, change in momt. assumed to be entirely due to forcing.
  //In this case, leave row diagonal to compute change in momt for post/dbg.
  //If dof (row) is free then i need to fill the non-diagonal terms.
  if( bForcedInSimFrame[0] ) { //then momenta not conserved in this dof
    A[0*6+1] = 0; A[0*6+2] = 0; A[0*6+3] = 0; A[0*6+4] = 0; A[0*6+5] = 0;
    b[0] = penalM * transVel_imposed[0]; // multply by penalM for conditioning
  }
  if( bForcedInSimFrame[1] ) { //then momenta not conserved in this dof
    A[1*6+0] = 0; A[1*6+2] = 0; A[1*6+3] = 0; A[1*6+4] = 0; A[1*6+5] = 0;
    b[1] = penalM * transVel_imposed[1];
  }
  if( bForcedInSimFrame[2] ) { //then momenta not conserved in this dof
    A[2*6+0] = 0; A[2*6+1] = 0; A[2*6+3] = 0; A[2*6+4] = 0; A[2*6+5] = 0;
    b[2] = penalM * transVel_imposed[2];
  }
  if( bBlockRotation[0] ) { //then momenta not conserved in this dof
    A[3*6+0] = 0; A[3*6+1] = 0; A[3*6+2] = 0; A[3*6+4] = 0; A[3*6+5] = 0;
    b[3] = 0; // TODO IMPOSED ANG VEL?
  }
  if( bBlockRotation[1] ) { //then momenta not conserved in this dof
    A[4*6+0] = 0; A[4*6+1] = 0; A[4*6+2] = 0; A[4*6+3] = 0; A[4*6+5] = 0;
    b[4] = 0; // TODO IMPOSED ANG VEL?
  }
  if( bBlockRotation[2] ) { //then momenta not conserved in this dof
    A[5*6+0] = 0; A[5*6+1] = 0; A[5*6+2] = 0; A[5*6+3] = 0; A[5*6+4] = 0;
    b[5] = 0; // TODO IMPOSED ANG VEL?
  }

  gsl_matrix_view Agsl = gsl_matrix_view_array (A.data(), 6, 6);
  gsl_vector_view bgsl = gsl_vector_view_array (b, 6);
  gsl_vector *xgsl = gsl_vector_alloc (6);
  int sgsl;
  gsl_permutation * permgsl = gsl_permutation_alloc (6);
  gsl_linalg_LU_decomp (& Agsl.matrix, permgsl, & sgsl);
  gsl_linalg_LU_solve (& Agsl.matrix, permgsl, & bgsl.vector, xgsl);
  transVel_computed[0] = gsl_vector_get(xgsl, 0);
  transVel_computed[1] = gsl_vector_get(xgsl, 1);
  transVel_computed[2] = gsl_vector_get(xgsl, 2);
  angVel_computed[0]   = gsl_vector_get(xgsl, 3);
  angVel_computed[1]   = gsl_vector_get(xgsl, 4);
  angVel_computed[2]   = gsl_vector_get(xgsl, 5);
  gsl_permutation_free (permgsl);
  gsl_vector_free (xgsl);

  force[0] = mass * (transVel_computed[0] - transVel[0]) / sim.dt;
  force[1] = mass * (transVel_computed[1] - transVel[1]) / sim.dt;
  force[2] = mass * (transVel_computed[2] - transVel[2]) / sim.dt;
  const std::array<Real,3> dAv = {
    (angVel_computed[0] - angVel[0]) / sim.dt,
    (angVel_computed[1] - angVel[1]) / sim.dt,
    (angVel_computed[2] - angVel[2]) / sim.dt
  };
  torque[0] = J[0] * dAv[0] + J[3] * dAv[1] + J[4] * dAv[2];
  torque[1] = J[3] * dAv[0] + J[1] * dAv[1] + J[5] * dAv[2];
  torque[2] = J[4] * dAv[0] + J[5] * dAv[1] + J[2] * dAv[2];

  if(bForcedInSimFrame[0]) {
    assert( std::fabs(transVel[0] - transVel_imposed[0]) < 1e-12 );
    transVel[0] = transVel_imposed[0];
  } else transVel[0] = transVel_computed[0];

  if(bForcedInSimFrame[1]) {
    assert( std::fabs(transVel[1] - transVel_imposed[1]) < 1e-12 );
    transVel[1] = transVel_imposed[1];
  } else transVel[1] = transVel_computed[1];

  if(bForcedInSimFrame[2]) {
    assert( std::fabs(transVel[2] - transVel_imposed[2]) < 1e-12 );
    transVel[2] = transVel_imposed[2];
  } else transVel[2] = transVel_computed[2];

  if( bBlockRotation[0] ) {
    assert( std::fabs(angVel[0] - 0) < 1e-12 );
    angVel[0] = 0;
  } else angVel[0] = angVel_computed[0];

  if( bBlockRotation[1] ) {
    assert( std::fabs(angVel[1] - 0) < 1e-12 );
    angVel[1] = 0;
  } else angVel[1] = angVel_computed[1];

  if( bBlockRotation[2] ) {
    assert( std::fabs(angVel[2] - 0) < 1e-12 );
    angVel[2] = 0;
  } else angVel[2] = angVel_computed[2];

  if (collision_counter > 0)
  {
      collision_counter -= sim.dt;
      transVel[0] = u_collision;
      transVel[1] = v_collision;
      transVel[2] = w_collision;
      angVel[0] = ox_collision;
      angVel[1] = oy_collision;
      angVel[2] = oz_collision;
  }

}

void Obstacle::computeForces()
{
  static const int nQoI = ObstacleBlock::nQoI;
  std::vector<Real> sum = std::vector<Real>(nQoI, 0);
  for (auto & block : obstacleBlocks) {
    if(block == nullptr) continue;
    block->sumQoI(sum);
  }

  MPI_Allreduce(MPI_IN_PLACE, sum.data(), nQoI, MPI_Real, MPI_SUM, sim.comm);

  //additive quantities: (check against order in sumQoI of ObstacleBlocks.h )
  unsigned k = 0;
  surfForce[0]  = sum[k++]; surfForce[1]  = sum[k++]; surfForce[2]  = sum[k++];
  presForce[0]  = sum[k++]; presForce[1]  = sum[k++]; presForce[2]  = sum[k++];
  viscForce[0]  = sum[k++]; viscForce[1]  = sum[k++]; viscForce[2]  = sum[k++];
  surfTorque[0] = sum[k++]; surfTorque[1] = sum[k++]; surfTorque[2] = sum[k++];
  drag          = sum[k++]; thrust        = sum[k++]; Pout          = sum[k++];
  PoutBnd       = sum[k++]; defPower      = sum[k++]; defPowerBnd   = sum[k++];
  pLocom        = sum[k++];

  const Real vel_norm = std::sqrt(transVel[0]*transVel[0]
                                  + transVel[1]*transVel[1]
                                  + transVel[2]*transVel[2]);
  //derived quantities:
  Pthrust    = thrust*vel_norm;
  Pdrag      =   drag*vel_norm;
  EffPDef    = Pthrust/(Pthrust-std::min(defPower,(Real)0)+EPS);
  EffPDefBnd = Pthrust/(Pthrust-         defPowerBnd        +EPS);

  _writeSurfForcesToFile();
  _writeDiagForcesToFile();
  _writeComputedVelToFile();
}

void Obstacle::update()
{
  const Real dqdt[4] = {
    (Real).5*( - angVel[0]*quaternion[1] - angVel[1]*quaternion[2] - angVel[2]*quaternion[3] ),
    (Real).5*( + angVel[0]*quaternion[0] + angVel[1]*quaternion[3] - angVel[2]*quaternion[2] ),
    (Real).5*( - angVel[0]*quaternion[3] + angVel[1]*quaternion[0] + angVel[2]*quaternion[1] ),
    (Real).5*( + angVel[0]*quaternion[2] - angVel[1]*quaternion[1] + angVel[2]*quaternion[0] )
  };

  if (sim.step < sim.step_2nd_start)
  {
    old_position  [0] = position  [0];
    old_position  [1] = position  [1];
    old_position  [2] = position  [2];
    old_absPos    [0] = absPos    [0];
    old_absPos    [1] = absPos    [1];
    old_absPos    [2] = absPos    [2];
    old_quaternion[0] = quaternion[0];
    old_quaternion[1] = quaternion[1];
    old_quaternion[2] = quaternion[2];
    old_quaternion[3] = quaternion[3];
    position  [0] += sim.dt * ( transVel[0] + sim.uinf[0] );
    position  [1] += sim.dt * ( transVel[1] + sim.uinf[1] );
    position  [2] += sim.dt * ( transVel[2] + sim.uinf[2] );
    absPos    [0] += sim.dt * transVel[0];
    absPos    [1] += sim.dt * transVel[1];
    absPos    [2] += sim.dt * transVel[2];
    quaternion[0] += sim.dt * dqdt[0];
    quaternion[1] += sim.dt * dqdt[1];
    quaternion[2] += sim.dt * dqdt[2];
    quaternion[3] += sim.dt * dqdt[3];
  }
  else
  {
    const Real aux = 1.0 / sim.coefU[0];

    Real temp [10] = {position[0],position[1],position[2],absPos[0],absPos[1],absPos[2],quaternion[0],quaternion[1],quaternion[2],quaternion[3]};
    position  [0] = aux * ( sim.dt * ( transVel[0] + sim.uinf[0] ) + ( - sim.coefU[1]*position  [0] - sim.coefU[2]*old_position  [0]) );
    position  [1] = aux * ( sim.dt * ( transVel[1] + sim.uinf[1] ) + ( - sim.coefU[1]*position  [1] - sim.coefU[2]*old_position  [1]) );
    position  [2] = aux * ( sim.dt * ( transVel[2] + sim.uinf[2] ) + ( - sim.coefU[1]*position  [2] - sim.coefU[2]*old_position  [2]) );
    absPos    [0] = aux * ( sim.dt * ( transVel[0]               ) + ( - sim.coefU[1]*absPos    [0] - sim.coefU[2]*old_absPos    [0]) );
    absPos    [1] = aux * ( sim.dt * ( transVel[1]               ) + ( - sim.coefU[1]*absPos    [1] - sim.coefU[2]*old_absPos    [1]) );
    absPos    [2] = aux * ( sim.dt * ( transVel[2]               ) + ( - sim.coefU[1]*absPos    [2] - sim.coefU[2]*old_absPos    [2]) );
    quaternion[0] = aux * ( sim.dt * ( dqdt[0]                   ) + ( - sim.coefU[1]*quaternion[0] - sim.coefU[2]*old_quaternion[0]) );
    quaternion[1] = aux * ( sim.dt * ( dqdt[1]                   ) + ( - sim.coefU[1]*quaternion[1] - sim.coefU[2]*old_quaternion[1]) );
    quaternion[2] = aux * ( sim.dt * ( dqdt[2]                   ) + ( - sim.coefU[1]*quaternion[2] - sim.coefU[2]*old_quaternion[2]) );
    quaternion[3] = aux * ( sim.dt * ( dqdt[3]                   ) + ( - sim.coefU[1]*quaternion[3] - sim.coefU[2]*old_quaternion[3]) );
    old_position  [0] = temp[0];
    old_position  [1] = temp[1];
    old_position  [2] = temp[2];
    old_absPos    [0] = temp[3];
    old_absPos    [1] = temp[4];
    old_absPos    [2] = temp[5];
    old_quaternion[0] = temp[6];
    old_quaternion[1] = temp[7];
    old_quaternion[2] = temp[8];
    old_quaternion[3] = temp[9];
  }
  const Real invD = 1.0/std::sqrt(quaternion[0]*quaternion[0] + quaternion[1]*quaternion[1] + quaternion[2]*quaternion[2] + quaternion[3]*quaternion[3]);
  quaternion[0] *= invD;
  quaternion[1] *= invD;
  quaternion[2] *= invD;
  quaternion[3] *= invD;


/*
  // normality preserving advection (Simulation of colliding constrained rigid bodies - Kleppmann 2007 Cambridge University, p51)
  // move the correct distance on the quaternion unit ball surface, end up with normalized quaternion
  const Real DQ[4] = { dqdt[0]*dt, dqdt[1]*dt, dqdt[2]*dt, dqdt[3]*dt };
  const Real DQn = std::sqrt(DQ[0]*DQ[0]+DQ[1]*DQ[1]+DQ[2]*DQ[2]+DQ[3]*DQ[3]);

  if(DQn>EPS)// && currentRKstep == 0)
  {
    const Real tanF = std::tan(DQn)/DQn;
    const Real D[4] = {
      Q[0] +tanF*DQ[0], Q[1] +tanF*DQ[1], Q[2] +tanF*DQ[2], Q[3] +tanF*DQ[3],
    };
    const Real invD = 1/std::sqrt(D[0]*D[0]+D[1]*D[1]+D[2]*D[2]+D[3]*D[3]);
    quaternion[0] = D[0] * invD; quaternion[1] = D[1] * invD;
    quaternion[2] = D[2] * invD; quaternion[3] = D[3] * invD;
  }
*/

  if (sim.verbose && sim.time > 0)
  {
    const Real rad2deg = 180./M_PI;
    std::array<Real,3> ypr = getYawPitchRoll();
    ypr[0] *= rad2deg;
    ypr[1] *= rad2deg;
    ypr[2] *= rad2deg;
    printf("pos:[%+.2f %+.2f %+.2f], u:[%+.2f %+.2f %+.2f], omega:[%+.2f %+.2f %+.2f], yaw: %+.1f, pitch: %+.1f, roll: %+.1f \n",
           absPos  [0],absPos  [1],absPos  [2],
           transVel[0],transVel[1],transVel[2],
           angVel  [0],angVel  [1],angVel  [2],
           ypr     [0],ypr     [1],ypr     [2]);
  }
  #ifndef NDEBUG
  const Real q_length=std::sqrt(quaternion[0]*quaternion[0]+quaternion[1]*quaternion[1]+quaternion[2]*quaternion[2]+quaternion[3]*quaternion[3]);
  assert(std::abs(q_length-1.0) < 5*EPS);
  #endif

}

void Obstacle::create()
{
  printf("Entered the wrong create operator\n");
  fflush(0); exit(1);
}

void Obstacle::finalize()
{ }

std::array<Real,3> Obstacle::getTranslationVelocity() const
{
  return std::array<Real,3> {{transVel[0],transVel[1],transVel[2]}};
}

std::array<Real,3> Obstacle::getAngularVelocity() const
{
  return std::array<Real,3> {{angVel[0],angVel[1],angVel[2]}};
}

std::array<Real,3> Obstacle::getCenterOfMass() const
{
  return std::array<Real,3> {{centerOfMass[0],centerOfMass[1],centerOfMass[2]}};
}

std::array<Real,3> Obstacle::getYawPitchRoll() const
{
  const Real roll  = atan2(      2.0 * (quaternion[3] * quaternion[2] + quaternion[0] * quaternion[1]) ,
                           1.0 - 2.0 * (quaternion[1] * quaternion[1] + quaternion[2] * quaternion[2]));
  const Real pitch = asin (      2.0 * (quaternion[2] * quaternion[0] - quaternion[3] * quaternion[1]));
  const Real yaw   = atan2(      2.0 * (quaternion[3] * quaternion[0] + quaternion[1] * quaternion[2]) ,
                          -1.0 + 2.0 * (quaternion[0] * quaternion[0] + quaternion[1] * quaternion[1]));
  return std::array<Real,3> {{yaw,pitch,roll}};
}

void Obstacle::saveRestart( FILE * f )
{
  assert(f != NULL);
  fprintf(f, "x:       %20.20e\n", (double)  position[0]);
  fprintf(f, "y:       %20.20e\n", (double)  position[1]);
  fprintf(f, "z:       %20.20e\n", (double)  position[2]);
  fprintf(f, "xAbs:    %20.20e\n", (double)    absPos[0]);
  fprintf(f, "yAbs:    %20.20e\n", (double)    absPos[1]);
  fprintf(f, "zAbs:    %20.20e\n", (double)    absPos[2]);
  fprintf(f, "quat_0:  %20.20e\n", (double)quaternion[0]);
  fprintf(f, "quat_1:  %20.20e\n", (double)quaternion[1]);
  fprintf(f, "quat_2:  %20.20e\n", (double)quaternion[2]);
  fprintf(f, "quat_3:  %20.20e\n", (double)quaternion[3]);
  fprintf(f, "u_x:     %20.20e\n", (double)  transVel[0]);
  fprintf(f, "u_y:     %20.20e\n", (double)  transVel[1]);
  fprintf(f, "u_z:     %20.20e\n", (double)  transVel[2]);
  fprintf(f, "omega_x: %20.20e\n", (double)    angVel[0]);
  fprintf(f, "omega_y: %20.20e\n", (double)    angVel[1]);
  fprintf(f, "omega_z: %20.20e\n", (double)    angVel[2]);
  fprintf(f, "old_position_0:    %20.20e\n", (double)old_position  [0]);
  fprintf(f, "old_position_1:    %20.20e\n", (double)old_position  [1]);
  fprintf(f, "old_position_2:    %20.20e\n", (double)old_position  [2]);
  fprintf(f, "old_absPos_0:      %20.20e\n", (double)old_absPos    [0]);
  fprintf(f, "old_absPos_1:      %20.20e\n", (double)old_absPos    [1]);
  fprintf(f, "old_absPos_2:      %20.20e\n", (double)old_absPos    [2]);
  fprintf(f, "old_quaternion_0:  %20.20e\n", (double)old_quaternion[0]);
  fprintf(f, "old_quaternion_1:  %20.20e\n", (double)old_quaternion[1]);
  fprintf(f, "old_quaternion_2:  %20.20e\n", (double)old_quaternion[2]);
  fprintf(f, "old_quaternion_3:  %20.20e\n", (double)old_quaternion[3]);
}

void Obstacle::loadRestart( FILE * f )
{
  assert(f != NULL);
  bool ret = true;
  double temp;
  ret = ret && 1==fscanf(f, "x:       %le\n", &temp);   position[0] = temp;
  ret = ret && 1==fscanf(f, "y:       %le\n", &temp);   position[1] = temp;
  ret = ret && 1==fscanf(f, "z:       %le\n", &temp);   position[2] = temp;
  ret = ret && 1==fscanf(f, "xAbs:    %le\n", &temp);     absPos[0] = temp;
  ret = ret && 1==fscanf(f, "yAbs:    %le\n", &temp);     absPos[1] = temp;
  ret = ret && 1==fscanf(f, "zAbs:    %le\n", &temp);     absPos[2] = temp;
  ret = ret && 1==fscanf(f, "quat_0:  %le\n", &temp); quaternion[0] = temp;
  ret = ret && 1==fscanf(f, "quat_1:  %le\n", &temp); quaternion[1] = temp;
  ret = ret && 1==fscanf(f, "quat_2:  %le\n", &temp); quaternion[2] = temp;
  ret = ret && 1==fscanf(f, "quat_3:  %le\n", &temp); quaternion[3] = temp;
  ret = ret && 1==fscanf(f, "u_x:     %le\n", &temp);   transVel[0] = temp;
  ret = ret && 1==fscanf(f, "u_y:     %le\n", &temp);   transVel[1] = temp;
  ret = ret && 1==fscanf(f, "u_z:     %le\n", &temp);   transVel[2] = temp;
  ret = ret && 1==fscanf(f, "omega_x: %le\n", &temp);     angVel[0] = temp;
  ret = ret && 1==fscanf(f, "omega_y: %le\n", &temp);     angVel[1] = temp;
  ret = ret && 1==fscanf(f, "omega_z: %le\n", &temp);     angVel[2] = temp;
  ret = ret && 1==fscanf(f, "old_position_0:    %le\n", &temp); old_position  [0] = temp;
  ret = ret && 1==fscanf(f, "old_position_1:    %le\n", &temp); old_position  [1] = temp;
  ret = ret && 1==fscanf(f, "old_position_2:    %le\n", &temp); old_position  [2] = temp;
  ret = ret && 1==fscanf(f, "old_absPos_0:      %le\n", &temp); old_absPos    [0] = temp;
  ret = ret && 1==fscanf(f, "old_absPos_1:      %le\n", &temp); old_absPos    [1] = temp;
  ret = ret && 1==fscanf(f, "old_absPos_2:      %le\n", &temp); old_absPos    [2] = temp;
  ret = ret && 1==fscanf(f, "old_quaternion_0:  %le\n", &temp); old_quaternion[0] = temp;
  ret = ret && 1==fscanf(f, "old_quaternion_1:  %le\n", &temp); old_quaternion[1] = temp;
  ret = ret && 1==fscanf(f, "old_quaternion_2:  %le\n", &temp); old_quaternion[2] = temp;
  ret = ret && 1==fscanf(f, "old_quaternion_3:  %le\n", &temp); old_quaternion[3] = temp;

  if( (not ret) )
  {
    printf("Error reading restart file. Aborting...\n");
    fflush(0); abort();
  }
  if (sim.verbose == 0 && sim.rank == 0)
    printf("Restarting Object.. x: %le, y: %le, z: %le, u_x: %le, u_y: %le, u_z: %le\n",(double)absPos[0],(double)absPos[1],(double)absPos[2],transVel[0],(double)transVel[1],(double)transVel[2]);
}

void Obstacle::_writeComputedVelToFile()
{
  if(sim.rank!=0 || sim.muteAll) return;
  std::stringstream ssR;
  ssR<<"velocity_"<<obstacleID<<".dat";
  std::stringstream &savestream = logger.get_stream(ssR.str());

  if(sim.step==0 && not printedHeaderVels) {
    printedHeaderVels = true;
    savestream<<"step time CMx CMy CMz quat_0 quat_1 quat_2 quat_3 vel_x vel_y vel_z angvel_x angvel_y angvel_z mass J0 J1 J2 J3 J4 J5 yaw pitch roll"<<std::endl;
  }

  const Real rad2deg = 180./M_PI;
  std::array<Real,3> ypr = getYawPitchRoll();
  ypr[0] *= rad2deg;
  ypr[1] *= rad2deg;
  ypr[2] *= rad2deg;

  savestream<<sim.step<<" ";
  savestream.setf(std::ios::scientific);
  savestream.precision(std::numeric_limits<float>::digits10 + 1);
  savestream <<sim.time
    <<" "<<absPos    [0]<<" "<<absPos    [1]<<" "<<absPos    [2]
    <<" "<<quaternion[0]<<" "<<quaternion[1]<<" "<<quaternion[2]<<" "<<quaternion[3]
    <<" "<<transVel  [0]<<" "<<transVel  [1]<<" "<<transVel  [2]
    <<" "<<angVel    [0]<<" "<<angVel    [1]<<" "<<angVel    [2]
    <<" "<<mass<<" "<<J[0]<<" "<<J[1]<<" "<<J[2]<<" "<<J[3]<<" "<<J[4]<<" "<<J[5]
    <<" "<<ypr[0]<<" "<<ypr[1]<<" "<<ypr[2]<<std::endl;
}

void Obstacle::_writeSurfForcesToFile()
{
  if(sim.rank!=0 || sim.muteAll) return;
  std::stringstream fnameF, fnameP;
  fnameF<<"forceValues_"<<obstacleID<<".dat";
  std::stringstream &ssF = logger.get_stream(fnameF.str());
  if(sim.step==0) {
    ssF<<"step time Fx Fy Fz torque_x torque_y torque_z FxPres FyPres FzPres FxVisc FyVisc FzVisc drag thrust"<<std::endl;
  }

  ssF << sim.step << " ";
  ssF.setf(std::ios::scientific);
  ssF.precision(std::numeric_limits<float>::digits10 + 1);
  ssF<<sim.time<<" "<<surfForce[0]<<" "<<surfForce[1]<<" "<<surfForce[2]
     <<" "<<surfTorque[0]<<" "<<surfTorque[1]<<" "<<surfTorque[2]<<" "
     <<presForce[0]<<" "<<presForce[1]<<" "<<presForce[2]<<" "<<viscForce[0]
     <<" "<<viscForce[1]<<" "<<viscForce[2]<<" "<<drag<<" "<<thrust<<std::endl;

  fnameP<<"powerValues_"<<obstacleID<<".dat";
  std::stringstream &ssP = logger.get_stream(fnameP.str());
  if(sim.step==0) {
    ssP<<"time Pthrust Pdrag PoutBnd Pout PoutNew defPowerBnd defPower EffPDefBnd EffPDef"<<std::endl;
  }
  ssP.setf(std::ios::scientific);
  ssP.precision(std::numeric_limits<float>::digits10 + 1);
  ssP<<sim.time<<" "<<Pthrust<<" "<<Pdrag<<" "<<PoutBnd<<" "<<Pout<<" "<<pLocom<<" "<<defPowerBnd<<" "<<defPower<<" "<<EffPDefBnd<<" "<<EffPDef<<std::endl;
}

void Obstacle::_writeDiagForcesToFile()
{
  if(sim.rank!=0 || sim.muteAll) return;
  std::stringstream fnameF;
  fnameF<<"forceValues_penalization_"<<obstacleID<<".dat";
  std::stringstream &ssF = logger.get_stream(fnameF.str());
  if(sim.step==0) {
    ssF << "step time mass force_x force_y force_z torque_x torque_y torque_z penalLmom_x penalLmom_y penalLmom_z penalAmom_x penalAmom_y penalAmom_z penalCM_x penalCM_y penalCM_z linVel_comp_x linVel_comp_y linVel_comp_z angVel_comp_x angVel_comp_y angVel_comp_z penalM penalJ0 penalJ1 penalJ2 penalJ3 penalJ4 penalJ5" << std::endl;
  }

  ssF << sim.step << " ";
  ssF.setf(std::ios::scientific);
  ssF.precision(std::numeric_limits<float>::digits10 + 1);
  ssF<<sim.time<<" "<<mass<<" "<<force[0]<<" "<<force[1]<<" "<<force[2]<<" "
     <<torque[0]<<" "<<torque[1]<<" "<<torque[2]
     <<" "<<penalLmom[0]<<" "<<penalLmom[1]<<" "<<penalLmom[2]
     <<" "<<penalAmom[0]<<" "<<penalAmom[1]<<" "<<penalAmom[2]
     <<" "<<penalCM[0]<<" "<<penalCM[1]<<" "<<penalCM[2]
     <<" "<<transVel_computed[0]<<" "<<transVel_computed[1]<<" "<<transVel_computed[2]
     <<" "<<angVel_computed[0]<<" "<<angVel_computed[1]<<" "<<angVel_computed[2]
     <<" "<<penalM<<" "<<penalJ[0]<<" "<<penalJ[1]<<" "<<penalJ[2]
     <<" "<<penalJ[3]<<" "<<penalJ[4]<<" "<<penalJ[5] <<std::endl;
}

CubismUP_3D_NAMESPACE_END
