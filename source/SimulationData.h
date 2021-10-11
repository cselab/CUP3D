//
//  CubismUP_3D
//  Copyright (c) 2018 CSE-Lab, ETH Zurich, Switzerland.
//  Distributed under the terms of the MIT license.
//
//  Created by Guido Novati (novatig@ethz.ch).
//

#ifndef CubismUP_3D_SimulationData_h
#define CubismUP_3D_SimulationData_h

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
class SpectralManip;

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

  // Step Counter, Time, and stopping criteria (0 means inactive)
  int step=0, nsteps=0;
  double time=0, endTime=0;

  // Current and Old Timestep
  double dt = 0;
  double dt_old = 0;

  // CFL number
  double CFL=0;

  // MPI
  MPI_Comm app_comm;
  int rank, nprocs;

  // Blocks Per Dimension XYZ
  int bpdx, bpdy, bpdz;

  // Start and Maximal Level of Refinement
  int levelStart, levelMax;

  // Refinement and Compression Tolerances
  double Rtol, Ctol;

  // Simulation Domain
  std::array<Real, 3> extent;  // Uniform grid by default.
  Real maxextent;

  // Resulting Maximal and Minimal Gridspacing
  Real hmin, hmax;

  // Velocity of Frame of Reference
  std::array<Real, 3> uinf = {{0, 0, 0}};

  // Flow Parameters
  double nu;

  // Penalisation Parameters
  double lambda, DLM=1;

  // Switch for Explicit/Implicit Penalisation
  bool bImplicitPenalization = true;

  // Initial conditions
  std::string initCond = "zero";
  std::string icFromH5 = "";

  // uMax Channel flow
  Real uMax_forced = 0;

  // Measured Umax
  Real uMax_measured = 0;

  // Time stepping
  // if step < step_2nd_start, explicit Euler steps are performed
  //(used to initialize u_{n-1} and u_n that are needed for 2nd order timestep)
  int step_2nd_start;
  double coefU[3] = {1.5,-2.0,0.5};
  int rampup;

  // Poisson solver
  double PoissonErrorTol;
  double PoissonErrorTolRel;

  // SGS
  std::string sgs = "";
  double cs = 0.0;
  double cs2mean = 0, cs2stdev = 0, nuSgsMean = 0, nuSgsStdev = 0;
  bool bComputeCs2Spectrum = false;
  double timeAnalysis = 0;

  // Indicator for collision
  bool bCollision = false;

  // Dump Settingas
  int freqDiagnostics = 0;
  bool bDump=false;
  bool verbose;
  bool muteAll;

  // output
  int statsFreq=1;
  int saveFreq=0;
  double saveTime=0, nextSaveTime=0;
  std::string path4serialization = "./";
  bool dumpP;
  bool dumpChi;
  bool dumpOmega,dumpOmegaX,dumpOmegaY,dumpOmegaZ;
  bool dumpVelocity,dumpVelocityX,dumpVelocityY,dumpVelocityZ;
  // flags assume value 0 for dirichlet/unbounded, 1 for periodic, 2 for wall
  BCflag BCx_flag = dirichlet, BCy_flag = dirichlet, BCz_flag = dirichlet;

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
#endif // CubismUP_3D_SimulationData_h
