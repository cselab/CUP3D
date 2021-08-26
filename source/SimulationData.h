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
  cubism::Profiler * profiler = nullptr;

  FluidGridMPI * grid = nullptr;
  inline std::vector<cubism::BlockInfo>& vInfo() {
    return grid->getBlocksInfo();
  }

  FluidGridMPIPoisson * gridPoisson = nullptr;
  inline std::vector<cubism::BlockInfo>& vInfoPoisson() {
    return gridPoisson->getBlocksInfo();
  }

  // Utility to update minimal and maximal gridspacing
  void updateH()
  {
    std::vector<cubism::BlockInfo> & myInfos = grid->getBlocksInfo();

    #pragma omp parallel for schedule(static) reduction(min : hmin,hmax)
    for(size_t i=0; i<myInfos.size(); i++)
    {
      const cubism::BlockInfo & info = myInfos[i];
      hmin = std::min(hmin,  info.h_gridpoint);
      hmax = std::min(hmax, -info.h_gridpoint);
    }
    Real res [2] = {hmin,hmax};
    MPI_Allreduce(MPI_IN_PLACE, & res, 2, MPI_REAL, MPI_MIN, app_comm);
    hmin  =  res[0];
    hmax  = -res[1];
  }

  AMR * amr;
  AMR2 * amr2;

  //The protagonist
  ObstacleVector * obstacle_vector = nullptr;
  //The antagonist
  std::vector<Operator*> pipeline;
  PoissonSolverAMR * pressureSolver = nullptr;
  SpectralManip * spectralManip = nullptr;
  // simulation status
  // nsteps==0 means that this stopping criteria is not active
  int step=0, nsteps=0;
  // endTime==0  means that this stopping criteria is not active
  double time=0, endTime=0;
  double dt = 0;

  // mpi
  MPI_Comm app_comm;
  int rank=-1, nprocs=-1;

  // grid
  int bpdx=-1, bpdy=-1, bpdz=-1;
  Real maxextent = 1;
  std::array<Real, 3> extent = {{1, 0, 0}};  // Uniform grid by default.
  bool bImplicitPenalization = true;
  Real hmin=0, hmax=0;
  int levelMax,levelStart;
  double Rtol,Ctol;

  // flow variables
  std::array<Real, 3> uinf = {{0, 0, 0}};
  double nu=0, CFL=0, lambda=-1, DLM=1;

  // initial conditions
  std::string initCond = "zero";
  std::string spectralIC = "";
  std::vector<double> initCondModes;
  std::vector<double> initCondSpectrum;
  std::string icFromH5 = "";
  double k0 = 0, tke0 = 0;

  // forcing
  bool bChannelFixedMassFlux = false;
  Real uMax_forced = 0, uMax_measured = 0;
  bool spectralForcing = false;
  double turbKinEn_target = 0; // read from settings
  double enInjectionRate = 0; // read from settings
  double actualInjectionRate = 0; // computed by specralManip, post processing

  // sgs
  std::string sgs = "";
  double cs = 0.0;

  // computed by SGS, for post processing:
  double cs2mean = 0, cs2stdev = 0, nuSgsMean = 0, nuSgsStdev = 0;
  bool bComputeCs2Spectrum = false;

  // analysis
  std::string analysis;
  double timeAnalysis = 0;
  int freqAnalysis = 0;
  double analysisTime=0, nextAnalysisTime=0;
  double grad_mean = 0, grad_std=0;

  // time stepping
  // if step < step_2nd_start, explicit Euler steps are performed
  //(used to initialize u_{n-1} and u_n that are needed for 2nd order timestep)
  int TimeOrder{1}; // =1 or =2
  int step_2nd_start{1};
  double coefU [3] = {1.5,-2.0,0.5};

  // analysis (channel)
  std::vector<Real> Ux_avg_tgt;
  std::vector<Real> kx_avg_tgt;
  std::vector<Real> Ux_avg_msr;
  std::vector<Real> kx_avg_msr;
  Real reTau = 0.0;

  // simulation settings
  int freqDiagnostics = 0;
  bool bDump=false;
  int rampup = 100;
  bool verbose=false;
  bool muteAll = false;
  double PoissonErrorTol = 1e-6;
  double PoissonErrorTolRel = 1e-4;

  // output
  int statsFreq=1;
  int saveFreq=0;
  double saveTime=0, nextSaveTime=0;
  std::string path4serialization = "./";
  std::string useSolver = "";
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
