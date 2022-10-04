//
//  Cubism3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#include "InitialConditions.h"
#include "ProcessHelpers.h"
#include "../Obstacles/ObstacleVector.h"
#include "../poisson/PoissonSolverBase.h"
#include <random>

class PoissonSolverBase;

CubismUP_3D_NAMESPACE_BEGIN
using namespace cubism;

namespace {

class KernelIC
{
 public:
  KernelIC(const Real u) {}
  void operator()(const BlockInfo& info, VectorBlock& block) const
  {
    for(int iz=0; iz<VectorBlock::sizeZ; ++iz)
    for(int iy=0; iy<VectorBlock::sizeY; ++iy)
    for(int ix=0; ix<VectorBlock::sizeX; ++ix) block(ix,iy,iz).clear();
  }
};

class KernelIC_taylorGreen
{
  const std::array<Real, 3> ext;
  const Real uMax;
  const Real a = 2*M_PI / ext[0], b = 2*M_PI / ext[1], c = 2*M_PI / ext[2];
  const Real A = uMax, B = - uMax * ext[1] / ext[0];
 public:
  KernelIC_taylorGreen(const std::array<Real, 3> &extents, const Real U): ext{extents}, uMax(U) {}
  void operator()(const BlockInfo& info, VectorBlock& block) const
  {
    for(int iz=0; iz<VectorBlock::sizeZ; ++iz)
    for(int iy=0; iy<VectorBlock::sizeY; ++iy)
    for(int ix=0; ix<VectorBlock::sizeX; ++ix)
    {
      block(ix,iy,iz).clear();
      Real p[3]; info.pos(p, ix, iy, iz);
      block(ix,iy,iz).u[0] = A*std::cos(a*p[0])*std::sin(b*p[1])*std::sin(c*p[2]);
      block(ix,iy,iz).u[1] = B*std::sin(a*p[0])*std::cos(b*p[1])*std::sin(c*p[2]);
    }
  }
};

class KernelIC_channel
{
  const int dir;
  const std::array<Real, 3> ext;
  const Real uMax, H = ext[dir], FAC = 4*uMax/H/H; // FAC = 0.5*G/mu
  //umax =  0.5*G/mu * 0.25*H*H
 public:
  KernelIC_channel(const std::array<Real, 3> &extents, const Real U, const int _dir):
    dir(_dir), ext{extents}, uMax(U) {}

  void operator()(const BlockInfo& info, VectorBlock& block) const
  {
    for(int iz=0; iz<VectorBlock::sizeZ; ++iz)
    for(int iy=0; iy<VectorBlock::sizeY; ++iy)
    for(int ix=0; ix<VectorBlock::sizeX; ++ix) {
      block(ix,iy,iz).clear();
      Real p[3]; info.pos(p, ix, iy, iz);
      block(ix,iy,iz).u[0] = FAC * p[dir] * (H-p[dir]);
    }
  }
};

class KernelIC_channelrandom
{
  const int dir;
  const std::array<Real, 3> ext;
  const Real uMax, H = ext[dir], FAC = 4*uMax/H/H; // FAC = 0.5*G/mu
  //const Real delta_tau = 5.0/180;
  //umax =  0.5*G/mu * 0.25*H*H
 public:
  KernelIC_channelrandom(const std::array<Real, 3> &extents, const Real U, const int _dir):
    dir(_dir), ext{extents}, uMax(U) {}

  void operator()(const BlockInfo& info, VectorBlock& block) const
  {
    std::random_device seed;
    std::mt19937 gen(seed());
    std::normal_distribution<Real> dist(0.0, 0.01);
    for(int iz=0; iz<VectorBlock::sizeZ; ++iz)
    for(int iy=0; iy<VectorBlock::sizeY; ++iy)
    for(int ix=0; ix<VectorBlock::sizeX; ++ix) {
      block(ix,iy,iz).clear();
      block(ix,iy,iz).u[0] = dist(gen);
    }
    //If we set block(ix,iy,iz).u = U*(1.0+dist(gen)) the compiler gives 
    //an annoying warning. Doing this slower approach with two loops makes
    //the warning disappear. This won't impact performance as it's done
    //onle once per simulation (initial conditions).
    for(int iz=0; iz<ScalarBlock::sizeZ; ++iz)
    for(int iy=0; iy<ScalarBlock::sizeY; ++iy)
    for(int ix=0; ix<ScalarBlock::sizeX; ++ix) {
      Real p[3]; info.pos(p, ix, iy, iz);
      const Real U = FAC * p[dir] * (H-p[dir]);
      block(ix,iy,iz).u[0] = U * ( block(ix,iy,iz).u[0] + 1.0 );
    }
  }
};

class IC_vorticity
{
 public:
  SimulationData& sim;
  const int Ncoil = 90;
  std::vector<Real> phi_coil;
  std::vector<Real> x_coil;
  std::vector<Real> y_coil;
  std::vector<Real> z_coil;

  std::vector<Real> dx_coil; //tangent vector 
  std::vector<Real> dy_coil; //tangent vector 
  std::vector<Real> dz_coil; //tangent vector 

  IC_vorticity(SimulationData & s) : sim(s)
  {
    phi_coil.resize(Ncoil);
    x_coil.resize(Ncoil);
    y_coil.resize(Ncoil);
    z_coil.resize(Ncoil);
    const int m = 2;
    const Real dphi = 2.0*M_PI / Ncoil;
    for (int i = 0; i < Ncoil; i++)
    {
      const Real phi = i * dphi;
      phi_coil[i] = phi;
      const Real R = 0.05*sin(m*phi);
      x_coil[i] = R* cos(phi)   + 1.0;
      y_coil[i] = R* sin(phi)   + 1.0;
      z_coil[i] = R* cos(m*phi) + 1.0;
    }

    dx_coil.resize(Ncoil);
    dy_coil.resize(Ncoil);
    dz_coil.resize(Ncoil);
    for (int i = 0; i < Ncoil; i++)
    {
      const Real phi = i * dphi;
      phi_coil[i] = phi;
      const Real R = 0.05*sin(m*phi);
      const Real dR = 0.05* m * cos(m*phi);
      const Real sinphi = sin(phi);
      const Real cosphi = cos(phi);
      dx_coil[i] = dR* cosphi - R * sinphi;
      dy_coil[i] = dR* sinphi + R * cosphi;
      dz_coil[i] = dR* cos(m*phi) - m * R * sin(m*phi);
      const Real norm = 1.0 / pow(dx_coil[i]*dx_coil[i]+
                                  dy_coil[i]*dy_coil[i]+
                                  dz_coil[i]*dz_coil[i]+1e-21,0.5);
      dx_coil[i] *= norm;
      dy_coil[i] *= norm;
      dz_coil[i] *= norm;
    }
  }

  int nearestCoil(const Real x, const Real y, const Real z)
  {
    int retval = -1;
    Real d = 1e10;
    for (int i = 0; i < Ncoil; i++)
    {
      const Real dtest = (x_coil[i] - x)*(x_coil[i] - x)+
                         (y_coil[i] - y)*(y_coil[i] - y)+
                         (z_coil[i] - z)*(z_coil[i] - z);
      if (dtest < d)
      {
        retval = i;
        d = dtest;
      }
    }
    return retval;
  }

  ~IC_vorticity() = default;

  void vort(const Real x, const Real y, const Real z, Real & omega_x, Real & omega_y, Real & omega_z)
  {
    const int idx = nearestCoil(x,y,z);
    const Real r2 = (x_coil[idx] - x)*(x_coil[idx] - x)+
                    (y_coil[idx] - y)*(y_coil[idx] - y)+
                    (z_coil[idx] - z)*(z_coil[idx] - z);
    const Real mag = 1.0/(r2+1)/(r2+1); 

    omega_x = mag * dx_coil[idx];
    omega_y = mag * dy_coil[idx];
    omega_z = mag * dz_coil[idx];
  }

  void run()
  {
    //1. Fill vel with vorticity values
    const int nz = VectorBlock::sizeZ;
    const int ny = VectorBlock::sizeY;
    const int nx = VectorBlock::sizeX;
    std::vector<BlockInfo>& velInfo = sim.velInfo();
    #pragma omp parallel for
    for(size_t i=0; i<velInfo.size(); i++)
    {
      Real p[3];
      VectorBlock & VEL = (*sim.vel)(i);
      for(int iz=0; iz<nz; ++iz)
      for(int iy=0; iy<ny; ++iy)
      for(int ix=0; ix<nx; ++ix)
      {
        velInfo[i].pos(p,ix,iy,iz);
        vort(p[0],p[1],p[2],VEL(ix,iy,iz).u[0],VEL(ix,iy,iz).u[1],VEL(ix,iy,iz).u[2]);
      }
    }

    //2. Compute curl(omega)
    //   Here we use the "ComputeVorticity" function from ProcessHelpers.h
    //   This computes the curl of whatever is stored in sim.vel (currently the vorticity field)
    //   and saves it to tmpV.
    {
      ComputeVorticity findOmega(sim);
      findOmega(0);
    }

    //3. Solve nabla^2 u = - curl(omega)
    std::shared_ptr<PoissonSolverBase> pressureSolver;
    pressureSolver = makePoissonSolver(sim);

    Real PoissonErrorTol    = sim.PoissonErrorTol;
    Real PoissonErrorTolRel = sim.PoissonErrorTolRel;
    sim.PoissonErrorTol    = 0; //we are solving this once, so we set tolerance = 0
    sim.PoissonErrorTolRel = 0; //we are solving this once, so we set tolerance = 0
    for (int d = 0 ; d < 3 ; d ++)
    {
      //3a. Fill RHS with -omega[d] and set initial guess to zero.
      #pragma omp parallel for
      for(size_t i=0; i<velInfo.size(); i++)
      {
        const VectorBlock & TMPV = (*sim.tmpV)(i);
        ScalarBlock & PRES = (*sim.pres)(i);
        ScalarBlock & LHS  = (*sim.lhs )(i);
        for(int iz=0; iz<nz; ++iz)
        for(int iy=0; iy<ny; ++iy)
        for(int ix=0; ix<nx; ++ix)
        {
          PRES(ix,iy,iz).s = 0.0;
          LHS (ix,iy,iz).s = -TMPV(ix,iy,iz).u[d];
        }
      }

      //3b. solve poisson equation
      pressureSolver->solve();

      //3c. Fill vel
      #pragma omp parallel for
      for(size_t i=0; i<velInfo.size(); i++)
      {
        VectorBlock & VEL = (*sim.vel)(i);
        ScalarBlock & PRES = (*sim.pres)(i);
        for(int iz=0; iz<nz; ++iz)
        for(int iy=0; iy<ny; ++iy)
        for(int ix=0; ix<nx; ++ix)
        {
          VEL(ix,iy,iz).u[d] = PRES(ix,iy,iz).s;
        }
      }
    }
    sim.PoissonErrorTol    = PoissonErrorTol   ; //recover tolerance for pressure projection
    sim.PoissonErrorTolRel = PoissonErrorTolRel; //recover tolerance for pressure projection
  }
};

}  // anonymous namespace

static void initialPenalization(SimulationData& sim, const Real dt)
{
  const std::vector<BlockInfo>& velInfo = sim.velInfo();
  for (const auto& obstacle : sim.obstacle_vector->getObstacleVector()) {
    using CHI_MAT = Real[CUP_BLOCK_SIZEZ][CUP_BLOCK_SIZEY][CUP_BLOCK_SIZEX];
    using UDEFMAT = Real[CUP_BLOCK_SIZEZ][CUP_BLOCK_SIZEY][CUP_BLOCK_SIZEX][3];
    // TODO: Refactor to use only one omp parallel.
    #pragma omp parallel
    {
      const auto& obstblocks = obstacle->getObstacleBlocks();
      const std::array<Real,3> centerOfMass = obstacle->getCenterOfMass();
      const std::array<Real,3> uBody = obstacle->getTranslationVelocity();
      const std::array<Real,3> omegaBody = obstacle->getAngularVelocity();

      #pragma omp for schedule(dynamic)
      for (size_t i = 0; i < velInfo.size(); ++i)
      {
        const BlockInfo& info = velInfo[i];
        const auto pos = obstblocks[info.blockID];
        if(pos == nullptr) continue;

        VectorBlock& b = (*sim.vel)(i);
        CHI_MAT & __restrict__ CHI = pos->chi;
        UDEFMAT & __restrict__ UDEF = pos->udef;

        for(int iz=0; iz<VectorBlock::sizeZ; ++iz)
        for(int iy=0; iy<VectorBlock::sizeY; ++iy)
        for(int ix=0; ix<VectorBlock::sizeX; ++ix)
        {
          Real p[3]; info.pos(p, ix, iy, iz);
          p[0]-=centerOfMass[0]; p[1]-=centerOfMass[1]; p[2]-=centerOfMass[2];
          const Real object_UR[3] = {
              (Real) omegaBody[1]*p[2] - (Real) omegaBody[2]*p[1],
              (Real) omegaBody[2]*p[0] - (Real) omegaBody[0]*p[2],
              (Real) omegaBody[0]*p[1] - (Real) omegaBody[1]*p[0]
          };
          const Real U_TOT[3] = {
              (Real)uBody[0] + object_UR[0] + UDEF[iz][iy][ix][0],
              (Real)uBody[1] + object_UR[1] + UDEF[iz][iy][ix][1],
              (Real)uBody[2] + object_UR[2] + UDEF[iz][iy][ix][2]
          };
          // what if multiple obstacles share a block??
          // let's plus equal and wake up during the night to stress about it
          b(ix,iy,iz).u[0] += CHI[iz][iy][ix] * ( U_TOT[0] - b(ix,iy,iz).u[0] );
          b(ix,iy,iz).u[1] += CHI[iz][iy][ix] * ( U_TOT[1] - b(ix,iy,iz).u[1] );
          b(ix,iy,iz).u[2] += CHI[iz][iy][ix] * ( U_TOT[2] - b(ix,iy,iz).u[2] );
        }
      }
    }
  }
}

void InitialConditions::operator()(const Real dt)
{
  if(sim.initCond == "zero") {
    if(sim.verbose) printf("[CUP3D] - Zero-values initial conditions.\n");
    run(KernelIC(0));
  }
  if(sim.initCond == "taylorGreen") {
    if(sim.verbose) printf("[CUP3D] - Taylor Green vortex initial conditions.\n");
    run(KernelIC_taylorGreen(sim.extents, sim.uMax_forced));
  }
  if(sim.initCond == "channelRandom")
  {
    if(sim.verbose) printf("[CUP3D] - Channel flow random initial conditions.\n");
    if( sim.BCx_flag == wall ) {
      printf("ERROR: channel flow must be periodic or dirichlet in x.\n");
      fflush(0); abort();
    }
    const bool channelY = sim.BCy_flag==wall, channelZ = sim.BCz_flag==wall;
    if( (channelY && channelZ) or (!channelY && !channelZ) ) {
      printf("ERROR: wrong channel flow BC in y or z.\n");
      fflush(0); abort();
    }
    const int dir = channelY? 1 : 2;
    run(KernelIC_channelrandom(sim.extents, sim.uMax_forced, dir));
  }
  if(sim.initCond == "channel")
  {
    if(sim.verbose) printf("[CUP3D] - Channel flow initial conditions.\n");
    if( sim.BCx_flag == wall ) {
      printf("ERROR: channel flow must be periodic or dirichlet in x.\n");
      fflush(0); abort();
    }
    const bool channelY = sim.BCy_flag==wall, channelZ = sim.BCz_flag==wall;
    if( (channelY && channelZ) or (!channelY && !channelZ) ) {
      printf("ERROR: wrong channel flow BC in y or z.\n");
      fflush(0); abort();
    }
    const int dir = channelY? 1 : 2;
    run(KernelIC_channel(sim.extents, sim.uMax_forced, dir));
  }
  if(sim.initCond == "vorticity") {
    if(sim.verbose) printf("[CUP3D] - Vorticity initial conditions.\n");
    IC_vorticity ic_vorticity(sim);
    ic_vorticity.run();
  }

  {
    std::vector<cubism::BlockInfo>& chiInfo  = sim.chiInfo();
    //zero fields, going to contain Udef:
    #pragma omp parallel for schedule(static)
    for(unsigned i=0; i<chiInfo.size(); i++)
    {
      ScalarBlock& CHI  = (*sim.chi )(i);
      ScalarBlock& PRES = (*sim.pres)(i);
      ScalarBlock& LHS  = (*sim.lhs )(i);
      VectorBlock& TMPV = (*sim.tmpV)(i);
      VectorBlock& VOLD = (*sim.vOld)(i);
      for(int iz=0; iz<ScalarBlock::sizeZ; ++iz)
      for(int iy=0; iy<ScalarBlock::sizeY; ++iy)
      for(int ix=0; ix<ScalarBlock::sizeX; ++ix) {
        CHI (ix,iy,iz).s = 0;
        PRES(ix,iy,iz).s = 0;
	      LHS (ix,iy,iz).s = 0;
	      TMPV(ix,iy,iz).u[0] = 0;
	      TMPV(ix,iy,iz).u[1] = 0;
	      TMPV(ix,iy,iz).u[2] = 0;
        VOLD(ix,iy,iz).u[0] = 0;
        VOLD(ix,iy,iz).u[1] = 0;
        VOLD(ix,iy,iz).u[2] = 0;
      }
    }
    //store deformation velocities onto tmp fields:
    initialPenalization(sim, dt);
  }
}

CubismUP_3D_NAMESPACE_END
