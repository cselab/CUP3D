//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#include "operators/ObstaclesUpdate.h"
#include "obstacles/ObstacleVector.h"
#include "utils/MatArrayMath.h"
#include <gsl/gsl_linalg.h>

CubismUP_3D_NAMESPACE_BEGIN
using namespace cubism;

using CHIMAT = Real[CUP_BLOCK_SIZE][CUP_BLOCK_SIZE][CUP_BLOCK_SIZE];
using UDEFMAT = Real[CUP_BLOCK_SIZE][CUP_BLOCK_SIZE][CUP_BLOCK_SIZE][3];
static constexpr Real EPS = std::numeric_limits<Real>::epsilon();
static constexpr Real DBLEPS = std::numeric_limits<double>::epsilon();

struct KernelIntegrateFluidMomenta : public ObstacleVisitor
{
  const double lambda, dt;
  ObstacleVector * const obstacle_vector;
  const cubism::BlockInfo * info_ptr = nullptr;
  inline double dvol(const cubism::BlockInfo&info, const int x, const int y, const int z) const {
    double h[3]; info.spacing(h, x, y, z);
    return h[0] * h[1] * h[2];
  }

  KernelIntegrateFluidMomenta(double _dt, double _lambda, ObstacleVector* ov)
    : lambda(_lambda), dt(_dt), obstacle_vector(ov) {}

  void operator()(const cubism::BlockInfo& info)
  {
    // first store the lab and info, then do visitor
    assert(info_ptr == nullptr);
    info_ptr = & info;
    ObstacleVisitor* const base = static_cast<ObstacleVisitor*> (this);
    assert( base not_eq nullptr );
    obstacle_vector->Accept( base );
    info_ptr = nullptr;
  }

  void visit(Obstacle* const op)
  {
    const BlockInfo& info = * info_ptr;
    assert(info_ptr not_eq nullptr);
    const std::vector<ObstacleBlock*>& obstblocks = op->getObstacleBlocks();
    ObstacleBlock*const o = obstblocks[info.blockID];
    if (o == nullptr) return;

    const std::array<double,3> CM = op->getCenterOfMass();
    const FluidBlock &b = *(FluidBlock *)info.ptrBlock;
    const CHIMAT & __restrict__ CHI = o->chi;
    const UDEFMAT & __restrict__ UDEF = o->udef;
    double &VV = o->V;
    double &FX = o->FX, &FY = o->FY, &FZ = o->FZ;
    double &TX = o->TX, &TY = o->TY, &TZ = o->TZ;
    VV = 0; FX = 0; FY = 0; FZ = 0; TX = 0; TY = 0; TZ = 0;
    double &J0 = o->J0, &J1 = o->J1, &J2 = o->J2;
    double &J3 = o->J3, &J4 = o->J4, &J5 = o->J5;
    J0 = 0; J1 = 0; J2 = 0; J3 = 0; J4 = 0; J5 = 0;
    o->GfX = 0;
    o->GpX = 0; o->GpY = 0; o->GpZ = 0;
    o->Gj0 = 0; o->Gj1 = 0; o->Gj2 = 0;
    o->Gj3 = 0; o->Gj4 = 0; o->Gj5 = 0;
    o->GuX = 0; o->GuY = 0; o->GuZ = 0;
    o->GaX = 0; o->GaY = 0; o->GaZ = 0;
    const Real lambdt = lambda*dt;
    for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
    for(int iy=0; iy<FluidBlock::sizeY; ++iy)
    for(int ix=0; ix<FluidBlock::sizeX; ++ix)
    {
      if (CHI[iz][iy][ix] <= 0) continue;
      double p[3]; info.pos(p, ix, iy, iz);
      const double dv = dvol(info, ix, iy, iz), X = CHI[iz][iy][ix];
      p[0] -= CM[0]; p[1] -= CM[1]; p[2] -= CM[2];

      VV += X * dv;
      J0 += X * dv * ( p[1]*p[1] + p[2]*p[2] );
      J1 += X * dv * ( p[0]*p[0] + p[2]*p[2] );
      J2 += X * dv * ( p[0]*p[0] + p[1]*p[1] );
      J3 -= X * dv * p[0]*p[1];
      J4 -= X * dv * p[0]*p[2];
      J5 -= X * dv * p[1]*p[2];

      FX += X * dv * b(ix,iy,iz).u;
      FY += X * dv * b(ix,iy,iz).v;
      FZ += X * dv * b(ix,iy,iz).w;
      TX += X * dv * ( p[1] * b(ix,iy,iz).w - p[2] * b(ix,iy,iz).v );
      TY += X * dv * ( p[2] * b(ix,iy,iz).u - p[0] * b(ix,iy,iz).w );
      TZ += X * dv * ( p[0] * b(ix,iy,iz).v - p[1] * b(ix,iy,iz).u );

      const Real penalFac = dv * lambdt * X / ( 1 + X * lambdt );
      o->GfX += penalFac;
      o->GpX += penalFac * p[0];
      o->GpY += penalFac * p[1];
      o->GpZ += penalFac * p[2];
      o->Gj0 += penalFac * ( p[1]*p[1] + p[2]*p[2] );
      o->Gj1 += penalFac * ( p[0]*p[0] + p[2]*p[2] );
      o->Gj2 += penalFac * ( p[0]*p[0] + p[1]*p[1] );
      o->Gj3 -= penalFac * p[0]*p[1];
      o->Gj4 -= penalFac * p[0]*p[2];
      o->Gj5 -= penalFac * p[1]*p[2];
      const double DiffU[3] = {
        b(ix,iy,iz).u - UDEF[iz][iy][ix][0],
        b(ix,iy,iz).v - UDEF[iz][iy][ix][1],
        b(ix,iy,iz).w - UDEF[iz][iy][ix][2]
      };
      o->GuX += penalFac * DiffU[0];
      o->GuY += penalFac * DiffU[1];
      o->GuZ += penalFac * DiffU[2];
      o->GaX += penalFac * ( p[1] * DiffU[2] - p[2] * DiffU[1] );
      o->GaY += penalFac * ( p[2] * DiffU[0] - p[0] * DiffU[2] );
      o->GaZ += penalFac * ( p[0] * DiffU[1] - p[1] * DiffU[0] );
    }
  }
};

struct KernelFinalizeObstacleVel : public ObstacleVisitor
{
  const double dt, lambda;
  FluidGridMPI * const grid;

  KernelFinalizeObstacleVel(double _dt, double _lambda, FluidGridMPI*g) :
    dt(_dt), lambda(_lambda), grid(g) { }

  void visit(Obstacle* const obst)
  {
    static constexpr int nQoI = 29;
    double M[nQoI] = { 0 };
    const auto& oBlock = obst->getObstacleBlocks();
    #pragma omp parallel for schedule(static,1) reduction(+ : M[:nQoI])
    for (size_t i=0; i<oBlock.size(); i++) {
      if(oBlock[i] == nullptr) continue;
      int k = 0;
      M[k++] += oBlock[i]->V ;
      M[k++] += oBlock[i]->FX; M[k++] += oBlock[i]->FY; M[k++] += oBlock[i]->FZ;
      M[k++] += oBlock[i]->TX; M[k++] += oBlock[i]->TY; M[k++] += oBlock[i]->TZ;
      M[k++] += oBlock[i]->J0; M[k++] += oBlock[i]->J1; M[k++] += oBlock[i]->J2;
      M[k++] += oBlock[i]->J3; M[k++] += oBlock[i]->J4; M[k++] += oBlock[i]->J5;
      M[k++] +=oBlock[i]->GfX;
      M[k++] +=oBlock[i]->GpX; M[k++] +=oBlock[i]->GpY; M[k++] +=oBlock[i]->GpZ;
      M[k++] +=oBlock[i]->Gj0; M[k++] +=oBlock[i]->Gj1; M[k++] +=oBlock[i]->Gj2;
      M[k++] +=oBlock[i]->Gj3; M[k++] +=oBlock[i]->Gj4; M[k++] +=oBlock[i]->Gj5;
      M[k++] +=oBlock[i]->GuX; M[k++] +=oBlock[i]->GuY; M[k++] +=oBlock[i]->GuZ;
      M[k++] +=oBlock[i]->GaX; M[k++] +=oBlock[i]->GaY; M[k++] +=oBlock[i]->GaZ;
      assert(k==29);
    }
    const auto comm = grid->getCartComm();
    MPI_Allreduce(MPI_IN_PLACE, M, nQoI, MPI_DOUBLE, MPI_SUM, comm);
    assert(std::fabs(obst->mass - M[ 0]) < 10*DBLEPS);
    assert(std::fabs(obst->J[0] - M[ 7]) < 10*DBLEPS);
    assert(std::fabs(obst->J[1] - M[ 8]) < 10*DBLEPS);
    assert(std::fabs(obst->J[2] - M[ 9]) < 10*DBLEPS);
    assert(std::fabs(obst->J[3] - M[10]) < 10*DBLEPS);
    assert(std::fabs(obst->J[4] - M[11]) < 10*DBLEPS);
    assert(std::fabs(obst->J[5] - M[12]) < 10*DBLEPS);
    assert(M[0] > DBLEPS);
    const double mass = M[0];
    // Obstacles may prevent rotation along certain coordinates, but we want to
    // still record the imposed angular acceleration for dbg/post.
    // Cross terms in the J matrix describe how angular momentum is transferred
    // across directions. If direction is locked it absorbs energy (ie. applies
    // constraint force) which prevents momentum transfer from that direction.
    // We accomplish this by zeroing some terms in invJ. E.g if obstacle can
    // only rotate around z, \dot{omega}_z = \tau_z / J_zz
    // If obstacle does not rotate in any direction we still want to see accels
    // (for debug) and therefore diagonal is never zeroed.
    const int bRotate[3] = { obst->bBlockRotation[0] ? 0 : 1,
                             obst->bBlockRotation[1] ? 0 : 1,
                             obst->bBlockRotation[2] ? 0 : 1 };
    const GenM fullJ = {
      M[ 7]             , M[10] * bRotate[1], M[11] * bRotate[2],
      M[10] * bRotate[0], M[ 8]             , M[12] * bRotate[2],
      M[11] * bRotate[0], M[12] * bRotate[1], M[ 9]
    };
    const GenM invJ = invertGen(fullJ);
    if(bRotate[0]==0 && bRotate[1]==0) // J_zz * invJzz == 1
      assert(std::fabs(invJ[8]*M[9]-1) < 10*DBLEPS);

    // these are the angular and linear momenta of the fluid inside the obstacle
    obst->angVel_fluid = multGenVec(invJ, GenV{{ M[4], M[5], M[6] }});
    obst->transVel_fluid = {M[1] / mass, M[2] / mass, M[3] / mass};

    #ifndef OLD_INTEGRATE_MOM
      const std::array<double,3> avobst = obst->getAngularVelocity();
      const std::array<double,3> utobst = obst->getTranslationVelocity();
      double A[6][6] = {0}, b[6] = {M[23], M[24], M[25], M[26], M[27], M[28]};
      //Momenta are conserved if a dof (a row of mat A) is not externally forced
      //This means that if obstacle is free to move according to fluid forces,
      //momenta after penal should be equal to moments before penal!
      //If dof is forced, change in momt. assumed to be entirely due to forcing.
      //In this case, leave row diagonal to compute change in momt for post/dbg.
      //If dof (row) is free then i need to fill the non-diagonal terms.
      if(obst->bForcedInSimFrame[0]) { //then momenta not conserved in this dof
        A[0][0] = M[13];
        b[0] = M[13] * utobst[0]; // multply by M_penal for conditioning
      } else {
        A[0][0] = M[13];
        A[0][1] = 0;
        A[0][2] = 0;
        A[0][3] = 0;
        A[0][4] = + M[16]; // PZ
        A[0][5] = - M[15]; // PY
      }
      if(obst->bForcedInSimFrame[1]) { //then momenta not conserved in this dof
        A[1][1] = M[13];
        b[1] = M[13] * utobst[1];
      } else {
        A[1][0] = 0;
        A[1][1] = M[13];
        A[1][2] = 0;
        A[1][3] = - M[16]; // PZ
        A[1][4] = 0;
        A[1][5] = + M[14]; // PX
      }
      if(obst->bForcedInSimFrame[2]) { //then momenta not conserved in this dof
        A[2][2] = M[13];
        b[2] = M[13] * utobst[2];
      } else {
        A[2][0] = 0;
        A[2][1] = 0;
        A[2][2] = M[13];
        A[2][3] = + M[15]; // PY
        A[2][4] = - M[14]; // PX
        A[2][5] = 0;
      }
      if(obst->bBlockRotation[0]) { //then momenta not conserved in this dof
        A[3][3] = M[17];
        b[3] = M[17] * avobst[0];
      } else {
        A[3][0] = 0;
        A[3][1] = - M[16]; // PZ
        A[3][2] = + M[15]; // PY
        A[3][3] = M[17]; // PJ0
        A[3][4] = M[20]; // PJ3
        A[3][5] = M[21]; // PJ4
      }
      if(obst->bBlockRotation[1]) { //then momenta not conserved in this dof
        A[4][4] = M[18];
        b[4] = M[18] * avobst[1];
      } else {
        A[4][0] = + M[16]; // PZ
        A[4][1] = 0;
        A[4][2] = - M[14]; // PX
        A[4][3] = M[20]; // PJ3
        A[4][4] = M[18]; // PJ1
        A[4][5] = M[22]; // PJ5
      }
      if(obst->bBlockRotation[2]) { //then momenta not conserved in this dof
        A[5][5] = M[19];
        b[5] = M[19] * avobst[2];
      } else {
        A[5][0] = - M[15]; // PY
        A[5][1] = + M[14]; // PX
        A[5][2] = 0;
        A[5][3] = M[21]; // PJ4
        A[5][4] = M[22]; // PJ5
        A[5][5] = M[19]; // PJ2
      }
      gsl_matrix_view Agsl = gsl_matrix_view_array (&A[0][0], 6, 6);
      gsl_vector_view bgsl = gsl_vector_view_array (b, 6);
      gsl_vector *xgsl = gsl_vector_alloc (6);
      int sgsl;
      gsl_permutation * permgsl = gsl_permutation_alloc (6);
      gsl_linalg_LU_decomp (& Agsl.matrix, permgsl, & sgsl);
      gsl_linalg_LU_solve (& Agsl.matrix, permgsl, & bgsl.vector, xgsl);
      obst->transVel_computed[0] = gsl_vector_get(xgsl, 0);
      obst->transVel_computed[1] = gsl_vector_get(xgsl, 1);
      obst->transVel_computed[2] = gsl_vector_get(xgsl, 2);
      obst->angVel_computed[0]   = gsl_vector_get(xgsl, 3);
      obst->angVel_computed[1]   = gsl_vector_get(xgsl, 4);
      obst->angVel_computed[2]   = gsl_vector_get(xgsl, 5);
      int rank;
      MPI_Comm_rank(MPI_COMM_WORLD,&rank);
      if(rank==0) {
        printf("Mom fluid: lin={%f %f %f}, ang={%f %f %f}\n Mom solid: lin={%f %f %f}, ang={%f %f %f}\n",
        obst->transVel_computed[0], obst->transVel_computed[1],
        obst->transVel_computed[2], obst->angVel_computed[0],
        obst->angVel_computed[1], obst->angVel_computed[2],
        obst->transVel_fluid[0], obst->transVel_fluid[1],
        obst->transVel_fluid[2], obst->angVel_fluid[0],
        obst->angVel_fluid[1], obst->angVel_fluid[2]); fflush(0);
      }
      obst->angVel_fluid = obst->angVel_computed;
      obst->transVel_fluid = obst->transVel_computed;
      gsl_permutation_free (permgsl);
      gsl_vector_free (xgsl);
      obst->force[0] = mass*(obst->transVel_computed[0]-obst->transVel[0])/dt;
      obst->force[1] = mass*(obst->transVel_computed[1]-obst->transVel[1])/dt;
      obst->force[2] = mass*(obst->transVel_computed[2]-obst->transVel[2])/dt;
      const std::array<double,3> dAv = {
        (obst->angVel_computed[0] - obst->angVel[0]) / dt,
        (obst->angVel_computed[1] - obst->angVel[1]) / dt,
        (obst->angVel_computed[2] - obst->angVel[2]) / dt
      };
      obst->torque[0] = M[ 7]*dAv[0] + M[10]*dAv[1] + M[11]*dAv[2];
      obst->torque[1] = M[10]*dAv[0] + M[ 8]*dAv[1] + M[12]*dAv[2];
      obst->torque[2] = M[11]*dAv[0] + M[12]*dAv[1] + M[ 9]*dAv[2];
    #else // momenta in the flow pre penal become momenta of the obstacle
      obst->force[0]  = M[1]; //these are all wrong now, must compute post penal
      obst->force[1]  = M[2]; //these are all wrong now, must compute post penal
      obst->force[2]  = M[3]; //these are all wrong now, must compute post penal
      obst->torque[0] = M[4]; //these are all wrong now, must compute post penal
      obst->torque[1] = M[5]; //these are all wrong now, must compute post penal
      obst->torque[2] = M[6]; //these are all wrong now, must compute post penal
      obst->transVel_computed[0] = obst->transVel_fluid[0];
      obst->transVel_computed[1] = obst->transVel_fluid[1];
      obst->transVel_computed[2] = obst->transVel_fluid[2];
      obst->angVel_computed[0] = obst->angVel_fluid[0];
      obst->angVel_computed[1] = obst->angVel_fluid[1];
      obst->angVel_computed[2] = obst->angVel_fluid[2];
    #endif

    obst->computeVelocities();
  }
};

void UpdateObstacles::operator()(const double dt)
{
  if(sim.obstacle_vector->nObstacles() == 0) return;

  sim.startProfiler("Obst Int Vel");
  { // integrate momenta by looping over grid
    #pragma omp parallel
    { // each thread needs to call its own non-const operator() function
      auto K = KernelIntegrateFluidMomenta(dt, sim.lambda, sim.obstacle_vector);
      #pragma omp for schedule(dynamic, 1)
      for (size_t i = 0; i < vInfo.size(); ++i) K(vInfo[i]);
    }
  }
  sim.stopProfiler();

  sim.startProfiler("Obst Upd Vel");
  {
    ObstacleVisitor*K = new KernelFinalizeObstacleVel(dt, sim.lambda, sim.grid);
    sim.obstacle_vector->Accept(K); // accept you son of a french cow
    delete K;
  }
  sim.stopProfiler();

  check("UpdateObstacles");
}

CubismUP_3D_NAMESPACE_END
