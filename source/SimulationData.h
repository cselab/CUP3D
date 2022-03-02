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
  FluidGridMPIPoisson * gridPoisson = nullptr;

  // Get velocity blocks on current rank
  inline std::vector<cubism::BlockInfo>& vInfo() {
    return grid->getBlocksInfo();
  }

  // Get pressure blocks on current rank
  inline std::vector<cubism::BlockInfo>& vInfoPoisson() {
    return gridPoisson->getBlocksInfo();
  }

  // Mesh Adaptation
  AMR * amr;
  AMR2 * amr2;

  // Container holding the obstacles
  ObstacleVector * obstacle_vector = nullptr;

  // Operator Pipeline
  std::vector<Operator*> pipeline;

  // Pressure solver to be shared between PressureRHS and PressureProjection
  PoissonSolverAMR * pressureSolver = nullptr;

  // Timestepping
  double dt      = 0;//current timestep
  double dt_old  = 0;//previous timestep
  double CFL     = 0;//Courant number
  double time    = 0;//current time
  int step       = 0;//currect step number
  double endTime = 0;//stop simulation at t=endTime (=0 means inactive)
  int nsteps     = 0;//stop simulation after nsteps (=0 means inactive)
  int rampup;        //exponential CFL rampup for the first 'rampup' steps
  int step_2nd_start;//explicit Euler for the first 'step_2nd_start' steps 
                     //(to initialize u_{n-1} for n=1)
  double coefU[3] = {1.5,-2.0,0.5};//used for 2nd order time integration 
                                   //of obstacle positions
  // MPI
  MPI_Comm app_comm;
  int rank, nprocs;

  //AMR & simulation domain
  int bpdx, bpdy, bpdz;                //blocks per dimension at refinement level 0
  int levelStart;                      //initial refinement level
  int levelMax;                        //max refinement level
  double Rtol;                         //mesh refinement tolerance
  double Ctol;                         //mesh compression tolerance
  std::array<double, 3> extent;        //simulation cubic domain extents
  double maxextent ;                   //max(extent[0],extent[1],extent[2])
  double hmin, hmax;                   //max and min grid spacing
  std::array<double, 3> uinf = {0,0,0};//velocity of Frame of Reference

  //Other stuff
  double uMax_measured = 0;         //max velocity at current timestep
  double nu;                        //fluid kinematic viscosity
  double lambda;                    //penalisation coefficient
  bool bImplicitPenalization = true;//explicit/implicit Penalisation
  double DLM=0;                     // if DLM>0 then lambda=DLM/dt
  double PoissonErrorTol;           //Poisson solver absolute error tolerance
  double PoissonErrorTolRel;        //Poisson solver relative error tolerance
  bool bCollision = false;          //indicator for collision between obstacles
  BCflag BCx_flag = freespace;      //boundary conditions in X
  BCflag BCy_flag = freespace;      //boundary conditions in Y
  BCflag BCz_flag = freespace;      //boundary conditions in Z

  // Initial conditions
  std::string initCond = "zero";
  std::string icFromH5 = "";

  // uMax Channel flow
  Real uMax_forced = 0;

  // Dump Settingas
  int freqDiagnostics = 0;
  int saveFreq = 0;
  bool bDump = false;
  bool verbose = false;
  bool muteAll = false;
  double saveTime=0;
  double nextSaveTime=0;
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
