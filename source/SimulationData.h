//
//  CubismUP_3D
//  Copyright (c) 2022 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//

#pragma once

#include "Definitions.h"
#include <array>
#include <vector>
#include <random>

namespace cubism {
  class Profiler;
  class ArgumentParser;
}

CubismUP_3D_NAMESPACE_BEGIN

class Operator;
class ObstacleVector;
class PoissonSolverAMR;

struct SimulationData
{
  // Profiler
  cubism::Profiler * profiler = nullptr;

  // Grids for Velocity and Pressure
  FluidGridMPI * grid = nullptr;
  ScalarGrid * lhs = nullptr;
  ScalarGrid * z   = nullptr;
  ScalarAMR * lhs_amr;
  ScalarAMR * z_amr;

  // Get velocity blocks on current rank
  inline std::vector<cubism::BlockInfo>& vInfo() {
    return grid->getBlocksInfo();
  }

  // Mesh Adaptation
  AMR * amr;

  // Container holding the obstacles
  ObstacleVector * obstacle_vector = nullptr;

  // Operator Pipeline
  std::vector<std::shared_ptr<Operator>> pipeline;

  // Pressure solver to be shared between PressureRHS and PressureProjection
  PoissonSolverAMR * pressureSolver = nullptr;

  // Timestepping
  Real dt      = 0;//current timestep
  Real dt_old  = 0;//previous timestep
  Real CFL     = 0;//Courant number
  Real time    = 0;//current time
  int step       = 0;//currect step number
  Real endTime = 0;//stop simulation at t=endTime (=0 means inactive)
  int nsteps     = 0;//stop simulation after nsteps (=0 means inactive)
  int rampup;        //exponential CFL rampup for the first 'rampup' steps
  int step_2nd_start;//explicit Euler for the first 'step_2nd_start' steps 
                     //(to initialize u_{n-1} for n=1)
  Real coefU[3] = {1.5,-2.0,0.5};//used for 2nd order time integration 
                                   //of obstacle positions
  // MPI
  MPI_Comm app_comm;
  int rank, nprocs;

  //AMR & simulation domain
  int bpdx, bpdy, bpdz;                //blocks per dimension at refinement level 0
  int levelStart;                      //initial refinement level
  int levelMax;                        //max refinement level
  Real Rtol;                         //mesh refinement tolerance
  Real Ctol;                         //mesh compression tolerance
  std::array<Real, 3> extent;        //simulation cubic domain extents
  Real maxextent ;                   //max(extent[0],extent[1],extent[2])
  Real hmin, hmax;                   //max and min grid spacing
  std::array<Real, 3> uinf = {0,0,0};//velocity of Frame of Reference

  //Other stuff
  Real uMax_measured = 0;         //max velocity at current timestep
  Real nu;                        //fluid kinematic viscosity
  Real lambda;                    //penalisation coefficient
  bool bImplicitPenalization = true;//explicit/implicit Penalisation
  Real DLM=0;                     // if DLM>0 then lambda=DLM/dt
  Real PoissonErrorTol;           //Poisson solver absolute error tolerance
  Real PoissonErrorTolRel;        //Poisson solver relative error tolerance
  bool bCollision = false;          //indicator for collision between obstacles
  BCflag BCx_flag = freespace;      //boundary conditions in X
  BCflag BCy_flag = freespace;      //boundary conditions in Y
  BCflag BCz_flag = freespace;      //boundary conditions in Z
  int bMeanConstraint = 1;          // if 0, no zero mean constraint for Poisson
                                    // if 1, replace one equation with zero mean contraint
				    // if 2, add mean to LHS. Ie. solve nabla^2 P + mean(P) = RHS
				    // which is mathematically equivalent to the zero-mean 
				    // solution of nabla^2 P = RHS (for zero Neumann BCs of course)
				    // if >2, we set one grid point to p=0

  // Initial conditions
  std::string initCond = "zero";
  std::string icFromH5 = "";

  // uMax Channel flow
  Real uMax_forced = 0;

  // Dump Settingas
  int freqDiagnostics = 0;
  int freqProfiler = 0;
  int saveFreq = 0;
  bool bDump = false;
  bool verbose = false;
  bool muteAll = false;
  Real saveTime=0;
  Real nextSaveTime=0;
  std::string path4serialization = "./";
  bool dumpP;
  bool dumpChi;
  bool dumpOmega,dumpOmegaX,dumpOmegaY,dumpOmegaZ;
  bool dumpVelocity,dumpVelocityX,dumpVelocityY,dumpVelocityZ;

  void startProfiler(std::string name) const;
  void stopProfiler() const;
  void printResetProfiler();
  void _preprocessArguments();
  ~SimulationData();
  SimulationData() = delete;
  SimulationData(const SimulationData &) = delete;
  SimulationData(SimulationData &&) = delete;
  SimulationData &operator=(const SimulationData &) = delete;
  SimulationData &operator=(SimulationData &&) = delete;
  SimulationData(MPI_Comm mpicomm, cubism::ArgumentParser &parser);
};

CubismUP_3D_NAMESPACE_END
